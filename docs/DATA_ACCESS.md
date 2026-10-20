# Data access

Competition files are listed by `GET https://api.virtualembryo.ai/challenge/data/manifest` and served through signed
links (`/challenge/data/link?key=...`, valid 15 minutes). Both need the bearer token the site keeps in the browser
after you sign in. Registration must be complete (Rules and Terms accepted, Human Team) or the manifest is refused.

## 1. Put the token on the local machine (never on the shared DGX account)

1. Sign in at https://virtualembryo.ai/challenge, open the browser DevTools console on that page and run:
   ```js
   copy(JSON.parse(localStorage.getItem("ve-challenge-session")).token)
   ```
   This copies the token to the clipboard.
2. In a terminal on the local machine:
   ```bash
   mkdir -p ~/.config/vec && chmod 700 ~/.config/vec
   cat > ~/.config/vec/token      # paste, press Enter, then Ctrl-D
   chmod 600 ~/.config/vec/token
   ```

The token expires after a while; if a command reports HTTP 401, repeat step 1.

## 2. Download onto the DGX

```bash
uv run vec data manifest            # list files and sizes
uv run vec data fetch               # every file -> PrakashDGX_H2:~/v/virtual-embryo-challenge/data/raw/
uv run vec data fetch --match T1    # a subset
```

`fetch` asks for one signed link at a time, right before that file, and hands it to `curl` on the DGX through stdin,
so neither the token nor the links appear in a process list or log. Complete files (size equal to the manifest's) are
skipped, so re-running after the final phase opens (2026-10-27) only pulls the new validation answers and test inputs.

Record SHA-256 checksums after each download (`sha256sum data/raw/**/*.h5ad > data/raw/SHA256SUMS` on the DGX) and
compare with the checksums the organisers publish on the Data page.
