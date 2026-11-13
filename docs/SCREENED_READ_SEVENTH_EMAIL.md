# Screened read for an answer to the seventh organiser email

Approved by the captain on 2026-11-13 for 19 and 23 Nov ("yes to the screened reads on 19 and 23 Nov"; `docs/PLAN.md`
decision log 2026-11-13 (d)); the check of 26-27 Nov repeats it. Written by the AI assistant before any read. The
main session never opens the sources below: it hands this file to a screening agent as its brief and reads only the
block that agent releases, after a second agent has audited that block against the scope.

## The brief (give it to the screening agent as written)

Do not open any file the injected memory index names.

You are a screening agent for a team in the Virtual Embryo Challenge. You read three sources for one thing only and
release nothing else. The team's main session must never see the rest.

Removed: test-board limits, anchors and truth counts, and any other property of a hidden target; per-board leaderboard
figures; gene-level knockout findings from outside the released data (the ruling of 2026-10-15 below). Kept: overall
totals and ranks; phenotype-level knowledge, which that ruling allows; our own results; the organisers' published
baselines.

**What the team asked the organisers on 2026-11-13** (subject "Human Team: processing batches of two held-out samples
(procedure only)"): (1) were the E7.75 sections (Task 2 embryo) processed together with those of any released MERFISH
stage, meaning on the same slide or hybridised and imaged in the same batch, and if so which stage; (2) were the E12.5
Multiome libraries (Task 1) prepared together with the released E10.5 library whose cell names end in "_3" and
sequenced in the same run, and if not, were they made and sequenced with the same protocol as that library.

**Sources.** Open these and nothing else:
1. The organisers' Data Protocol page, https://virtualembryo.ai/challenge/data-protocol. Not the Data page, not the
   boards, not the leaderboard, no file linked from it. If it answers HTTP 403 or asks for a sign-in, stop and report
   that; no other user agent, no other address, no cached copy.
2. The public Slack channel #virtual-embryo-challenge: organiser posts since 2026-11-13 only. Skip any message that
   quotes test-board numbers, veckit 0.1.2 or PR #6.
3. The organiser thread in the captain's mailbox: replies to the subject above, and nothing else in the mailbox.

**Release only this:**
- For each source: "no answer found", or the organisers' statement, quoted, that says which released sample or stage
  the E7.75 sections or the E12.5 Multiome libraries were processed with, or that they were processed separately, or
  that a protocol was the same or another. Replace every number in a quote with [n], stage names (E7.5, E8.0, E10.5,
  E12.5 and the like) and the library suffix "_3" excepted.
- For the Data Protocol page: whether its text about how the hidden E7.75 and E12.5 samples were processed has changed
  since the read of 2026-11-10, yes or no, without the changed text unless it is a statement of the kind above.
- The date and time of each read.

**Withhold everything else,** and say only how many passages you withheld and of which class: a measured or relative
property of any held-out sample (a size, a count, a depth, a share, a number of sections or cells, a file size);
anything about a held-out sample other than the two asked about; test-board limits or formats; scores; gene-level
statements about a knockout; anything you are unsure of.

**Your final answer** is the released block only, under 200 words, plus the withheld counts by class.

## What the main session does with the released block

Record it with its date in `docs/ORGANIZER_QUESTIONS.md` ("Seventh email") and `docs/DATA_POLICY.md`, and apply the
reading fixed there beforehand. Nothing in an upload changes without the captain's ruling. If a reply reached the
captain's inbox first, the captain says so and does not paste it; this read then covers it.
