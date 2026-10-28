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
skipped, so re-running only pulls what is new. No test inputs are distributed.

In the final phase the P2 validation answers (T1 E10.5, T2 E8.5, E10.5 and E7.5, the two T3 Gata4 replicates) have
their own list and links (`/challenge/data/validation/state`, `.../manifest`, `.../link?key=...`), served only after
the team has unlocked them. A person unlocks them on the site's Data page, after ticking "I understand and agree to
permanently stop P2 submissions for my entire team"; it cannot be undone, even without a download, and no tool in
this repo sends it. Our team unlocked them on 2026-10-28 at 09:32 UTC. `vec data manifest` prints the unlock state
and both lists, and `vec data fetch` downloads both into the same `data/raw/T*/` tree (the keys do not overlap; a key
listed in both is refused).

Record SHA-256 checksums after each download (`sha256sum data/raw/**/*.h5ad > data/raw/SHA256SUMS` on the DGX) and
compare sizes with the manifest: the site publishes no checksums (as of 2026-10-28). The 15 files on the DGX are
listed in `docs/DATA_NOTES.md`.
