# Mathematical derivations — mega27-08-phage-design

## 1. Pair-AUROC for additive scoring rules (corrected statement)

**Setup.** Each phage $i$ appears twice: once with its true host panel
$(x_i, c^{+}_i, y{=}1)$ and once with the alternative panel
$(x_i, c^{-}_i, y{=}0)$. Let $A$ be the phages whose true host is species 1
(panel $c_1$), $B$ those whose true host is species 2 (panel $c_2$),
$n_A = |A|$, $n_B = |B|$.

**Theorem (balanced case).** If the scoring rule is additive,
$s(x, c) = f(x) + g(c)$, and $n_A = n_B$, then the pair-AUROC equals $1/2$
for every $f, g$.

**Proof.** Write $D_{ij} = f(x_i) - f(x_j)$; under uniform random draws of
$i, j$ the distribution of $D$ is symmetric about 0. Partition the
positive-negative comparison pairs by host class of $i$ and $j$:

* $i \in A, j \in B$ or $i \in B, j \in A$ (cross-class): the host terms
  cancel ($g(c_1)$ vs $g(c_1)$, or $g(c_2)$ vs $g(c_2)$), so each comparison
  is decided by $D_{ij} > 0$, contributing exactly $1/2$.
* $i, j \in A$: decided by $D_{ij} > g(c_2) - g(c_1) =: -\delta$,
  contributing $P(D > -\delta)$.
* $i, j \in B$: decided by $D_{ij} > +\delta$, contributing
  $P(D > \delta) = 1 - P(D > -\delta)$ (ties have measure zero).

The same-class cells contribute
$\tfrac{n_A^2 P(D>-\delta) + n_B^2 (1 - P(D>-\delta))}{(n_A+n_B)^2}$ and
the cross-class cells $\tfrac{2 n_A n_B \cdot 1/2}{(n_A+n_B)^2}$. When
$n_A = n_B$ the same-class term is also $1/2$, so AUROC $= 1/2$. $\blacksquare$

**Corollary (imbalance caveat).** With $n_A \neq n_B$ an additive rule can
score above or below $1/2$ (by up to
$\tfrac{|n_A^2 - n_B^2|}{(n_A+n_B)^2} \cdot |P(D>-\delta) - \tfrac12|$),
so the degeneracy statement must not be quoted as exact for imbalanced
panels; our dataset (431 vs 359 phages) is mildly imbalanced, and the
qualitative point - additive rules cannot express phage-specific host
preference - is what excludes them as baselines.

**Retracted empirical claim.** An earlier draft cited an observed baseline
AUROC of exactly 0.500 as empirical confirmation. That run's fit had in fact
collapsed (lbfgs stopped at the initial point, $n_{iter}=0$, coefficient norm
0), so the observation was an optimization artifact, not evidence for the
theorem. It is preserved here as a caution: *verify convergence before
reading anything into a metric.* After fixing feature scaling, the additive
baseline was dropped on the Corollary's grounds, and only coupling models
(Hadamard k-mer, bilinear, hybrid) are benchmarked.

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
