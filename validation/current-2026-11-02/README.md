# Actual current-publication validation

The exact declared runtime closure was copied away from the development and
publication working trees. All mathematical processes used isolated Python
flags (-I, -S, -B), explicitly allowed installed libraries, file-read/output
guards and disabled network access. The four mathematical child scripts were
run unchanged under the same guards. No private source or publication checkout
was read, and all251 copied source files remained unchanged.

The [actual public-copy run](public-copy-replay.json) read all163 definitions and
all163 integral entries, then ran the complete calculation:

| Stage | Actual checks | Record |
| --- | --- | --- |
| Low-loop individual-I color/scales/CT/RG | 51, including2 negative controls | [low](low-checks.json) |
| Three matter colors from31 native orbits | 101, including2 negative controls | [matter](matter-checks.json) |
| AA from current source/joint objects | 30, including2 negative controls, plus217 absolute-order convolutions | [AA](AA-checks.json) |
| Combined four-color scales/CT/RG and exact function representatives | 197, including2 negative controls | [combined](combine-checks.json) |

The newly generated one/two/three-loop result files were byte-identical to the
published current result files. These counts include related algebraic checks;
they are not counts of independent physical proofs.

The [earlier independent reader review](original-reader-review.json) additionally
records492 checks including17 refusal controls on the same API/profile bytes:
mixed orbit/normalization contracts, extra signs, duplicate Euler compensation,
wrong orbit/loop references, changed hashes and escaping paths. Batch reader
coverage caches a single validated immutable snapshot; negative controls restore
the full loader. The fresh publication run above separately executes the full
mathematical replay, including the final AA output guard.

The [current manifest](../../workflow/current-calculation.json) binds runtime
source/code/data. The separate [publication manifest](../../CURRENT-PUBLICATION.json)
binds this dated export and its English documentation. Older upstream and low-loop
manifests remain tied to their own snapshots; mutable current navigation is not
silently rebound to those full-tree histories.

Original high-dimensional integrations, function dictionaries and admitted
research solutions are explicitly reused. Hashes, successful software calls and
RG alone do not prove the original physics or fix an arbitrary scale-independent
finite shape. No new CAS integration, arbitrary complex-domain theorem or Lean
whole-integral proof is claimed. formal=not_run.
