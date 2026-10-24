# Security Policy

## Supported versions

| Version | Supported |
|---------|-----------|
| 1.0.x   | ✅ Current — deployed at plantacrispai.netlify.app |
| < 1.0   | ❌ Pre-release, unsupported |

## The threat model

PlantaCRISPAI is a **fully client-side application**. There is no backend, no
database, no authentication and no server-side session. Sequence data is parsed,
aligned and classified in the user's browser and is never transmitted.

This removes most classes of vulnerability, and concentrates the remaining risk
in three places:

1. **The `.ab1` binary parser.** It consumes untrusted binary input from an
   arbitrary file. Malformed or hostile input causing memory exhaustion, an
   infinite loop, or a crash is the most plausible defect class here.
2. **Rendering of file-derived strings.** Sample names and run dates come from
   the uploaded file and are displayed in the UI.
3. **The supply chain.** A compromised npm dependency reaching the deployed
   bundle would be the highest-impact failure, since it would run in the context
   of a page handling unpublished sequence data.

### What would be a vulnerability

- Any path by which sequence data, file contents or file names leave the client
- Script execution from a crafted `.ab1`, FASTA or pasted sequence
- A crafted input that hangs or crashes the tab (denial of service)
- Dependency vulnerabilities reaching the production bundle
- Missing or bypassable security headers on the deployment

### What would not be

- Incorrect genotype calls or scientific inaccuracy — important, but a
  [bug report](https://github.com/shelvys/plantacrispai/issues), not a security
  issue
- Results differing from ICE or TIDE — see
  [`docs/SCIENTIFIC_METHODS.md`](docs/SCIENTIFIC_METHODS.md)
- Findings against `react-scripts` dev-only dependencies that are not in the
  production bundle
- Missing rate limiting, authentication or account controls — there is no server
  and no account

## Reporting a vulnerability

**Please do not open a public issue for a security problem.**

Report privately to **<crispricar@gmail.com>** with subject line
`SECURITY: PlantaCRISPAI`.

Please include:

- a description of the issue and why you believe it is exploitable
- steps to reproduce, with a minimal proof of concept if you have one
- affected version or deployment URL
- browser and operating system
- any suggested remediation

> [!CAUTION]
> If a proof of concept requires a `.ab1` file, send a synthetic or already
> published one. Do not send unpublished research data.

### What to expect

| Stage | Target |
|---|---|
| Acknowledgement of your report | within 5 working days |
| Initial assessment and severity | within 10 working days |
| Fix for a confirmed high-severity issue | within 30 days |
| Public disclosure | coordinated with you, after a fix is deployed |

We will credit you in the release notes and in
[`CHANGELOG.md`](CHANGELOG.md) unless you prefer to remain anonymous.

Please allow a reasonable period for remediation before public disclosure. We
have no bug bounty programme — this is a publicly funded academic research tool —
but genuine findings are gratefully received and credited.

## Deployment hardening

The production deployment sets, via [`netlify.toml`](netlify.toml):

- `X-Content-Type-Options: nosniff`
- `X-Frame-Options: SAMEORIGIN`
- `Referrer-Policy: strict-origin-when-cross-origin`
- Immutable, long-lived caching restricted to content-hashed `/static/*` assets

HTTPS is enforced by the hosting platform.

## Guidance for users handling sensitive data

- Sequence data is processed locally, but your **browser** is still a trust
  boundary. Keep it updated and be mindful of extensions with page access.
- The app works offline once loaded. For maximum assurance with unpublished
  data, load the page, disconnect from the network, then analyse.
- Downloaded reports and CSV exports are ordinary files on your machine and are
  no longer protected by the app's privacy model.
- Browser history and the page title do not record sequence content.

## Dependencies

Runtime dependencies are deliberately limited to React and React DOM, which
minimises supply-chain surface. Dependency updates are monitored by Dependabot;
see [`.github/dependabot.yml`](.github/dependabot.yml).
