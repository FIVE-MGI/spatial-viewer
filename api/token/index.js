// POST /api/token  {password}  -> {base, sas, expires}
// GET  /api/token             -> same, but authorised by the session cookie
//
// The SpatialViewer data lives in a PRIVATE blob container, so a URL on its own
// fetches nothing. This endpoint is the only thing that knows the storage key:
// it checks the shared password and returns a container-scoped, read-only SAS
// that expires after TOKEN_HOURS. The key itself never reaches the browser.
//
// App settings (Azure portal -> Static Web App -> Environment variables):
//   VIEWER_PASSWORD    the shared password people type
//   STORAGE_KEY        access key for the storage account; mints a fresh SAS per
//                      session. Preferred.
//   BLOB_SAS           a ready-made read-only SAS, used only when STORAGE_KEY is
//                      absent - for when we have Blob Data Contributor (which cannot
//                      read keys) rather than key access.
//   STORAGE_ACCOUNT    required
//   STORAGE_CONTAINER  required
//   TOKEN_HOURS        default 1
//   SESSION_HOURS      default 1    (how long a browser can renew without retyping)
const crypto = require("crypto");
const {
  StorageSharedKeyCredential,
  generateBlobSASQueryParameters,
  ContainerSASPermissions,
  SASProtocol,
} = require("@azure/storage-blob");

const ACCOUNT = process.env.STORAGE_ACCOUNT || "";
const CONTAINER = process.env.STORAGE_CONTAINER || "";
// An hour each. The SAS is deliberately not longer-lived than the session that
// obtained it: a token copied out of devtools should not outlive the right to have
// asked for it, and three days of silent renewal made "change the password" a far
// weaker action than it looks.
const TOKEN_HOURS = Number(process.env.TOKEN_HOURS || 1);
const SESSION_HOURS = Number(process.env.SESSION_HOURS || 1);
const COOKIE = "sv_session";

// Brute force throttle. In-memory, so it is per running instance and resets on
// restart - enough to make guessing slow, not a substitute for a good password.
// Deliberately generous: a whole lab behind one campus IP or VPN looks like a
// single address here, and a tripped lockout rejects the right password too, so
// a handful of typos must not lock everyone else out. 30 tries / 5 min is still
// only ~8,600 a day, which gets a guesser nowhere.
const MAX_FAILS = 30;
const LOCKOUT_MS = 5 * 60 * 1000;
const fails = new Map();

function clientIp(req) {
  const fwd = req.headers["x-forwarded-for"] || "";
  return String(fwd).split(",")[0].trim() || "unknown";
}
function locked(ip) {
  const f = fails.get(ip);
  if (!f) return false;
  if (Date.now() > f.until) { fails.delete(ip); return false; }
  return f.count >= MAX_FAILS;
}
function noteFail(ip) {
  const f = fails.get(ip) || { count: 0, until: 0 };
  f.count += 1;
  f.until = Date.now() + LOCKOUT_MS;
  fails.set(ip, f);
}
function sameString(a, b) {
  const x = Buffer.from(String(a), "utf8");
  const y = Buffer.from(String(b), "utf8");
  if (x.length !== y.length) return false;         // length alone is not secret here
  return crypto.timingSafeEqual(x, y);
}

// ---- session cookie: "<expiry ms>.<hmac>", signed with the storage key -------
// It carries no data beyond its own expiry; it only says "this browser typed the
// password at some point", which is what lets a long-open tab renew its SAS
// without asking again.
function secret() {
  // Signed with the password as well as the server-side key. Without the password in
  // here, changing VIEWER_PASSWORD revoked nothing: every browser already holding a
  // cookie went on renewing its SAS for the rest of SESSION_HOURS, so the one action
  // an administrator takes to cut off access did not cut off anybody. Folding the
  // password in means changing it invalidates every outstanding session, which is
  // what changing a password is for.
  const key = process.env.STORAGE_KEY || process.env.BLOB_SAS || "";
  const pw = process.env.VIEWER_PASSWORD || "";
  return crypto.createHash("sha256").update("sv-session|" + key + "|" + pw).digest();
}
function signSession(expiresAt) {
  const mac = crypto.createHmac("sha256", secret()).update(String(expiresAt)).digest("base64url");
  return `${expiresAt}.${mac}`;
}
function validSession(value) {
  if (!value) return false;
  const dot = String(value).lastIndexOf(".");
  if (dot < 1) return false;
  const expiresAt = String(value).slice(0, dot);
  if (!/^\d+$/.test(expiresAt) || Number(expiresAt) < Date.now()) return false;
  return sameString(String(value), signSession(expiresAt));
}
function readCookie(req, name) {
  const raw = req.headers.cookie || "";
  for (const part of raw.split(";")) {
    const eq = part.indexOf("=");
    if (eq > 0 && part.slice(0, eq).trim() === name) return part.slice(eq + 1).trim();
  }
  return null;
}

function mintSas() {
  // Fallback for when we hold no account key. Storage Blob Data Contributor can move
  // data but deliberately cannot read account keys, so a SAS made elsewhere (a user
  // delegation SAS, which Azure caps at 7 days) can be pasted into BLOB_SAS instead.
  // The password still gates everything; what is lost is the per-session 8 hour
  // expiry, because every visitor now gets the same token until that SAS lapses.
  const preset = (process.env.BLOB_SAS || "").replace(/^\?/, "").trim();
  if (!process.env.STORAGE_KEY && preset) {
    const se = new URLSearchParams(preset).get("se");
    const expires = se ? Date.parse(se) : Date.now() + TOKEN_HOURS * 3600 * 1000;
    return {
      base: `https://${ACCOUNT}.blob.core.windows.net/${CONTAINER}/`,
      sas: preset,
      expires: Number.isFinite(expires) ? expires : Date.now() + 3600 * 1000,
    };
  }

  const key = process.env.STORAGE_KEY;
  if (!key) throw new Error("neither STORAGE_KEY nor BLOB_SAS is set");
  const cred = new StorageSharedKeyCredential(ACCOUNT, key);
  const expiresOn = new Date(Date.now() + TOKEN_HOURS * 3600 * 1000);
  const sas = generateBlobSASQueryParameters(
    {
      containerName: CONTAINER,
      permissions: ContainerSASPermissions.parse("r"),   // read only; no list, no write, no delete
      startsOn: new Date(Date.now() - 5 * 60 * 1000),    // 5 min back for clock skew
      expiresOn,
      protocol: SASProtocol.Https,
    },
    cred
  ).toString();
  return {
    base: `https://${ACCOUNT}.blob.core.windows.net/${CONTAINER}/`,
    sas,
    expires: expiresOn.getTime(),
  };
}

module.exports = async function (context, req) {
  const noStore = { "Content-Type": "application/json", "Cache-Control": "no-store" };
  const ip = clientIp(req);

  // A tab that already has a session just gets a fresh SAS.
  if (req.method === "GET") {
    if (!validSession(readCookie(req, COOKIE))) {
      context.res = { status: 401, headers: noStore, body: { error: "no session" } };
      return;
    }
    try {
      context.res = { status: 200, headers: noStore, body: mintSas() };
    } catch (e) {
      context.log.error("SAS failed:", e.message);
      context.res = { status: 500, headers: noStore, body: { error: "server not configured" } };
    }
    return;
  }

  if (locked(ip)) {
    context.res = { status: 429, headers: noStore, body: { error: "too many attempts - wait a few minutes" } };
    return;
  }

  const expected = process.env.VIEWER_PASSWORD;
  if (!expected) {
    context.log.error("VIEWER_PASSWORD app setting is not set");
    context.res = { status: 500, headers: noStore, body: { error: "server not configured" } };
    return;
  }

  const body = req.body || {};
  const given = typeof body === "string" ? (() => { try { return JSON.parse(body).password; } catch (e) { return ""; } })() : body.password;

  if (!given || !sameString(given, expected)) {
    noteFail(ip);
    await new Promise((r) => setTimeout(r, 400));     // slow down guessing
    context.res = { status: 401, headers: noStore, body: { error: "wrong password" } };
    return;
  }
  fails.delete(ip);

  let out;
  try {
    out = mintSas();
  } catch (e) {
    context.log.error("SAS failed:", e.message);
    context.res = { status: 500, headers: noStore, body: { error: "server not configured" } };
    return;
  }

  const sessionUntil = Date.now() + SESSION_HOURS * 3600 * 1000;
  const cookie = [
    `${COOKIE}=${signSession(sessionUntil)}`,
    "HttpOnly",
    "Secure",
    "SameSite=Strict",
    "Path=/api",
    `Max-Age=${SESSION_HOURS * 3600}`,
  ].join("; ");

  context.res = { status: 200, headers: { ...noStore, "Set-Cookie": cookie }, body: out };
};
