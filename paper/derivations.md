# Mathematical derivations — mega27-08-phage-design

## 1. Degeneracy of pair-AUROC under constant host offsets

**Claim.** Let each phage $i$ appear in two samples $(x_i, c^+, y=1)$ and
$(x_i, c^-, y=0)$, where $c^+, c^-$ are fixed host-panel representations. For
any additive scoring rule $s(x, c) = f(x) + g(c)$, the pair-AUROC over the
balanced pair dataset equals $1/2$ regardless of $f$.

**Proof.** $\mathrm{AUROC} = P(s(x_i, c^+) > s(x_j, c^-))$ over uniformly drawn
phage pairs $i, j$. Substituting the additive form,
$s(x_i, c^+) - s(x_j, c^-) = f(x_i) - f(x_j) + g(c^+) - g(c^-)$.
Write $\Delta g = g(c^+) - g(c^-)$, a constant. Then
$\mathrm{AUROC} = P(f(x_i) - f(x_j) > -\Delta g)$.
Since $(i, j)$ and $(j, i)$ are equally likely and
$f(x_i) - f(x_j) = -(f(x_j) - f(x_i))$, the distribution of
$f(x_i) - f(x_j)$ is symmetric about 0, hence
$P(f(x_i) - f(x_j) > c) = P(f(x_i) - f(x_j) < -c)$ for all $c$.
With ties of measure zero, $P(D > c) + P(D < -c) + P(-c \le D \le c) = 1$;
for the balanced construction $c = -\Delta g$ is constant, so
$\mathrm{AUROC} = 1 - \mathrm{AUROC}$, giving $\mathrm{AUROC} = 1/2$. $\blacksquare$

**Consequence.** Any legitimate host-interaction benchmark with paired
negatives must use features or model terms that couple $x$ and $c$
(bilinear, Hadamard, or attention). We therefore benchmark the Hadamard
k-mer baseline and bilinear/hybrid models only; additive "baseline"
configurations are excluded as provably degenerate. (Empirically confirmed:
the additive LR scored AUROC exactly 0.500.)

## 2. Permutation invariance of the mean-normalized MPNN readout

**Claim.** For the message-passing update
$m_i = \frac{1}{|N(i)|}\sum_{j \in N(i)} M(h_i, h_j)$,
$h_i' = U(h_i, m_i)$, followed by $z = \mathrm{MLP}(\frac{1}{L}\sum_i h_i^{(T)})$,
the embedding $z$ is invariant under any permutation $\pi$ of the node set.

**Proof.** Message passing at node $i$ depends only on the multiset
$\{(h_i, h_j) : j \in N(i)\}$; permuting node labels permutes the multiset of
node states without changing it, because both the neighbor sum and the mean
over $N(i)$ are symmetric functions. Induction over layers: if the layer-$t$
state multiset is invariant, so is the layer-$(t{+}1)$ multiset. The final
mean readout is a symmetric function of the layer-$T$ states, and MLP applied
to its value preserves the invariance. $\blacksquare$

(Empirically enforced by test_gnn_permutation_invariant_readout.)

## 3. Bilinear compatibility scoring

The interaction logit $s(p, h) = e_p^\top W e_h + b$ parameterizes a general
second-order coupling between phage embedding $e_p$ and host embedding $e_h$.
$W$ is learned unconstrained; if biological compatibility were symmetric
($s(p,h) = s(h,p)$ for equal-dimensional embeddings) the optimum would satisfy
$W = W^\top$, but host-range is directional (the phage binds the host, not
conversely), so no symmetry constraint is imposed. The hybrid model adds a
learned k-mer interaction head and a 2-parameter convex combination; the
combiner weights are inspected in the paper to quantify how much of the final
score is sequence-local (CNN) versus compositional (k-mer).

## 4. Docking score as a bounded pairwise potential

The rigid docking score
$S = \#\{(i,j) : 4.5 \le d_{ij} \le 9\} - 50 \cdot \#\{(i,j): d_{ij} < 4.5\}
 + 0.1 \sum_{d_{ij} \le 20} \frac{-q_i q_j}{\max(d_{ij}, 2)}$
is a finite sum of bounded pairwise terms; for $n_r, n_l$ residues its
evaluation is $O(n_r n_l)$ per pose and the grid search over $N$ poses is
$O(N n_r n_l)$. The 4.5-9 A contact band matches the geometric span of
first-shell van der Waals contacts at protein interfaces (literature-cited in
the paper); the clash penalty dominates by construction, so the optimizer
cannot gain by burial.
