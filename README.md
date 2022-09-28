# delta-classification

Sage code and corresponding classification data implementing the algorithms described in the paper [On the maximal number of columns of a \Delta-modular matrix](https://arxiv.org/abs/2111.06294). This is joint work with Gennadiy Averkov.

## Content:

### Code:

* **polytopes.sage**:               basic operations on lattice polytopes
* **delta-classification.sage**:    classification of o-symmetric lattice polytopes of a given dimension
                                    and a given value corresponding to the volume of a largest inscribed crosspolytope

### Data:

This folder includes lists of o-symmetric lattice polytopes P with a given dimension m and a given number \Delta = \Delta(P) (see Section 5 of the paper for the definition). Either full enumeration data is available for the pair (m,\Delta), or the enumeration data concerns all extremal examples, that is, those polytopes that maximize the number of integer points contained in P.

* **dim_2_delta_[1-13].txt**: o-symmetric lattice polytopes with \Delta(P) in [1-13] in dimension 2
* **dim_2_delta_[1-13]_extremal.txt**: extremal o-symmetric lattice polytopes with \Delta(P) in [1-13] in dimension 2
* **dim_3_delta_[1-7].txt**: o-symmetric lattice polytopes with \Delta(P) in [1-7] in dimension 3
* **dim_3_delta_[1-12]_extremal.txt**: extremal o-symmetric lattice polytopes with \Delta(P) in [1-12] in dimension 3
* **dim_4_delta_[1-2].txt**: o-symmetric lattice polytopes with \Delta(P) in [1-2] in dimension 4
* **dim_4_delta_[1-3]_extremal**: extremal o-symmetric lattice polytopes with \Delta(P) in [1-3] in dimension 4
