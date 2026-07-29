---
layout: post
title: "Nevanlinna-Pick Interpolation and RKHS"
date: 2026-07-29
---

# Nevanlinna-Pick interpolation

The classical Nevanlinna-Pick theorem is often presented as a matrix criterion. Given distinct points 

$$ z_{1},\ldots,z_{n}\in\mathbb{D} $$

and prescribed values 

$$ w_{1},\ldots,w_{n}\in\mathbb{D} $$

one asks whether there is a holomorphic function 

$$ f\colon\mathbb{D}\longrightarrow\mathbb{D} $$

such that 

$$ f\left(z_{i}\right)=w_{i}. $$

The theorem says that 

$$ f\;\text{exists}\Longleftrightarrow\left[\frac{1-w_{i}\overline{w_{j}}}{1-z_{i}\overline{z_{j}}}\right]^{n}_{i,j=1}\geq0. $$

The reproducing kernel Hilbert space (RKHS) formulation below gives a direct answer to why this holds. The interpolation values prescribe the action of the adjoint of a multiplier on finitely many kernel sections. The Pick matrix is the matrix of the condition that this prescribed action be contractive.

---

## The Szegő kernel and multiplier eigenvectors

Let $H^{2}$ be the Hardy space on the unit disk. Its reproducing kernel is the Szegő kernel 

$$ k\left(z,w\right)=\frac{1}{1-z\overline{w}}. $$

We write kernels sections as $k_{w}=k\left(\cdot,w\right)$, so that 

$$ h\left(w\right)=\left\langle k_{w},h\right\rangle _{H^{2}} $$

for every $h\in H^{2}$.

Suppose that $f$ is a bounded holomorphic function on $\mathbb{D}$, and let $M_{f}$ be multiplication by $f$ on $H^{2}$. The kernel functions are eigenvectors for the adjoint multiplier. Indeed, 

$$ M^{*}_{f}k_{w}=\overline{f\left(w\right)}k_{w}. $$

To see this, let $h\in H^{2}$. Then 

$$ \begin{aligned}\left\langle M^{*}_{f}k_{w},h\right\rangle  & =\left\langle k_{w},M_{f}h\right\rangle =f\left(w\right)h\left(w\right)=\left\langle \overline{f\left(w\right)}k_{w},h\right\rangle .\end{aligned} $$

Now suppose that we want $f\left(z_{i}\right)=w_{i}$. Any interpolating function must satisfy 

$$ M^{*}_{f}k_{z_{i}}=\overline{w_{i}}k_{z_{i}}. $$

The interpolation data therefore prescribe the restriction of $M^{*}_{f}$ to the finite-dimensional space 

$$ \mathcal{M}=\text{span}\left\{ k_{z_{1}},\ldots,k_{z_{n}}\right\} \subset H^{2}. $$

Since the points $z_{1},\ldots,z_{n}$ are distinct, the corresponding kernel functions are linearly independent. One way to see this is to use the power series expansion 

$$ k_{z_{i}}\left(z\right)=\frac{1}{1-\overline{z_{i}}z}=\sum^{\infty}_{m=0}\overline{z_{i}}^{m}z^{m}. $$

Suppose that $\sum^{n}_{i=1}c_{i}k_{z_{i}}=0$. Then, as an identity of holomorphic functions on $\mathbb{D}$, 

$$ \begin{aligned}0 & =\sum^{n}_{i=1}c_{i}k_{z_{i}}\left(z\right)=\sum^{n}_{i=1}c_{i}\sum^{\infty}_{m=0}\overline{z_{i}}^{m}z^{m}=\sum^{\infty}_{m=0}\left(\sum^{n}_{i=1}c_{i}\overline{z_{i}}^{m}\right)z^{m}.\end{aligned} $$

So, for all $m\geq0$, $\sum^{n}_{i=1}c_{i}\overline{z_{i}}^{m}=0$. Taking $m=0,\ldots,n-1$, this gives 

$$ \begin{pmatrix}1 & 1 & \cdots & 1\\
\overline{z_{1}} & \overline{z_{2}} & \cdots & \overline{z_{n}}\\
\vdots & \vdots &  & \vdots\\
\overline{z_{1}}^{n-1} & \overline{z_{2}}^{n-1} & \cdots & \overline{z_{n}}^{n-1}
\end{pmatrix}\begin{pmatrix}c_{1}\\
c_{2}\\
\vdots\\
c_{n}
\end{pmatrix}=0. $$

It's a Vandermonde matrix. Since the points $z_{1},\ldots,z_{n}$ are distinct, its determinant is 

$$ \prod_{1\leq i<j\leq n}\left(\overline{z_{j}}-\overline{z_{i}}\right)\neq0. $$

It follows that $c_{1}=\cdots=c_{n}=0$. Thus the kernel functions $k_{z_{1}},\ldots,k_{z_{n}}$ are linearly independent. 

We may therefore define an operator $A\colon\mathcal{M}\longrightarrow\mathcal{M}$ by 

$$ Ak_{z_{i}}=\overline{w_{i}}k_{z_{i}}. $$

The interpolation problem can now be written as an operator extension problem. We ask whether $A$ is the restriction of $M^{*}_{f}$ for some contractive multiplier $M_{f}$ on $H^{2}$.

---

## The Pick matrix as a contractivity condition

If $f$ maps $\mathbb{D}$ into $\overline{\mathbb{D}}$, then $\left\Vert M_{f}\right\Vert =\left\Vert f\right\Vert _{\infty}\leq1$. So $M^{*}_{f}$ is contractive, and its restriction $A$ must also be contractive.

Let $h=\sum^{n}_{i=1}c_{i}k_{z_{i}}\in\mathcal{M}$. Using $\left\langle k_{z_{i}},k_{z_{j}}\right\rangle =\frac{1}{1-z_{i}\overline{z_{j}}}$, we get 

$$ \left\Vert h\right\Vert ^{2}=\sum^{n}_{i,j=1}\overline{c_{i}}c_{j}\frac{1}{1-z_{i}\overline{z_{j}}}. $$

On the other hand, $Ah=\sum^{n}_{i=1}c_{i}\overline{w_{i}}k_{z_{i}}$, and so 

$$ \left\Vert Ah\right\Vert ^{2}=\sum^{n}_{i,j=1}\overline{c_{i}}c_{j}\frac{w_{i}\overline{w_{j}}}{1-z_{i}\overline{z_{j}}}. $$

It follows that 

$$ \left\Vert h\right\Vert ^{2}-\left\Vert Ah\right\Vert ^{2}=\sum^{n}_{i,j=1}\overline{c_{i}}c_{j}\frac{1-w_{i}\overline{w_{j}}}{1-z_{i}\overline{z_{j}}}. $$

Therefore 

$$ A\text{ is contractive}\Longleftrightarrow\left[\frac{1-w_{i}\overline{w_{j}}}{1-z_{i}\overline{z_{j}}}\right]^{n}_{i,j=1}\geq0. $$

This is the Pick matrix, i.e., the matrix of the quadratic form $I_{\mathcal{M}}-A^{*}A$. Its positivity says that the action prescribed on the kernel vectors is contractive. This gives one useful formulation of the theorem.

> **Theorem (Nevanlinna-Pick).**  
> Let $z_{1},\ldots,z_{n}$ be distinct points of $\mathbb{D}$, and let $w_{1},\ldots,w_{n}\in\mathbb{D}$. The following conditions are equivalent.
> 1. There is a holomorphic function $f\colon\mathbb{D}\longrightarrow\mathbb{D}$ such that $f\left(z_{i}\right)=w_{i}$.
> 2. The operator on $\mathcal{M}=\text{span}\left\{ k_{z_{1}},\ldots,k_{z_{n}}\right\}$ defined by $Ak_{z_{i}}=\overline{w_{i}}k_{z_{i}}$ is contractive.
> 3. The Pick matrix is positive semidefinite, i.e., 
> $$ \left[\frac{1-w_{i}\overline{w_{j}}}{1-z_{i}\overline{z_{j}}}\right]^{n}_{i,j=1}\geq0. $$

*Proof.* The equivalence between (2) and (3) is the calculation above. An interpolating Schur function gives (2) because $A$ is a restriction of $M^{*}_{f}$.

The remaining assertion is the extension theorem. For the Szegő kernel, every contractive action of this form on a finite kernel span comes from a contractive multiplier on the whole Hardy space. $\blacksquare$

---

## The de Branges-Rovnyak kernel

There is a second way to read the same matrix.

For a holomorphic function $f\colon\mathbb{D}\to\mathbb{D}$, define 

$$ c_{f}\left(z,w\right)=\frac{1-f\left(z\right)\overline{f\left(w\right)}}{1-z\overline{w}}. $$

This is the de Branges-Rovnyak kernel associated with $f$.

Its positivity follows directly from the multiplier inequality. Since $M_{f}$ is contractive, $I-M_{f}M^{*}_{f}\geq0$. Moreover, 

$$ \begin{aligned}\left\langle k_{z},\left(I-M_{f}M^{*}_{f}\right)k_{w}\right\rangle  & =\left\langle k_{z},k_{w}\right\rangle -\left\langle M^{*}_{f}k_{z},M^{*}_{f}k_{w}\right\rangle \\ & =\frac{1-f\left(z\right)\overline{f\left(w\right)}}{1-z\overline{w}}. \end{aligned} $$

Thus $c_{f}\left(z,w\right)=\left\langle k_{z},\left(I-M_{f}M^{*}_{f}\right)k_{w}\right\rangle$. The RKHS associated with $c_{f}$ is denoted by $\mathcal{H}\left(f\right)$.

If $f\left(z_{i}\right)=w_{i}$, then 

$$ c_{f}\left(z_{i},z_{j}\right)=\frac{1-w_{i}\overline{w_{j}}}{1-z_{i}\overline{z_{j}}}. $$

The Pick matrix is then the restriction of the unknown de Branges-Rovnyak kernel to the interpolation nodes.

Let $F=\left\{ z_{1},\ldots,z_{n}\right\}$ and define a kernel on $F$ by 

$$ \Gamma\left(z_{i},z_{j}\right)=\frac{1-w_{i}\overline{w_{j}}}{1-z_{i}\overline{z_{j}}}. $$

The Nevanlinna-Pick theorem can now be restated as follows.

> **Theorem.**  
> The interpolation data admit a Schur-class interpolant if and only if the kernel $\Gamma$ is positive definite.
> 
> Equivalently, $\Gamma$ is the restriction to $F\times F$ of a de Branges-Rovnyak kernel $c_{f}$ for some Schur function $f$.

This should not be confused with the general problem of extending a positive definite kernel from a finite set. The required extension must have the special form $\frac{1-f\left(z\right)\overline{f\left(w\right)}}{1-z\overline{w}}$ for one holomorphic function $f$. The theorem says that positivity of the finite kernel is sufficient for such an extension.

Thus the two RKHS formulations are complementary. The operator formulation says that the interpolation values prescribe a contractive action on a finite span of Szegő kernels. The kernel formulation says that the Pick matrix is the proposed finite Gram matrix of the de Branges-Rovnyak space of the interpolant.

---

## The lurking isometry

The sufficient part of the Pick condition can be proved using finite-dimensional Hilbert space geometry.

Assume that the Pick matrix is positive semidefinite. Choose a finite-dimensional Hilbert space $E$ and vectors $v_{1},\ldots,v_{n}\in E$ such that 

$$ \left\langle v_{i},v_{j}\right\rangle _{E}=\frac{1-w_{j}\overline{w_{i}}}{1-z_{j}\overline{z_{i}}}. $$

The order of the indices is chosen to agree with the convention that inner products are linear in the second variable.

The Gram identity can be rearranged as 

$$ 1+\overline{z_{i}}z_{j}\left\langle v_{i},v_{j}\right\rangle =\overline{w_{i}}w_{j}+\left\langle v_{i},v_{j}\right\rangle . $$

It follows that the two families of vectors 

$$ \begin{pmatrix}1\\ z_{i}v_{i} \end{pmatrix}\qquad\text{and}\qquad\begin{pmatrix}w_{i}\\ v_{i} \end{pmatrix} $$

in $\mathbb{C}\oplus E$ have the same Gram matrix.

Therefore the correspondence 

$$ \begin{pmatrix}1\\ z_{i}v_{i} \end{pmatrix}\longmapsto\begin{pmatrix}w_{i}\\ v_{i} \end{pmatrix} $$

extends linearly to an isometry between their spans. Since $\mathbb{C}\oplus E$ is finite dimensional, this isometry can be extended to a unitary operator $U\colon\mathbb{C}\oplus E\longrightarrow\mathbb{C}\oplus E$.

Write $U$ in block form as 

$$ U=\begin{pmatrix}A & B\\ C & D \end{pmatrix}. $$

Here $A\in\mathbb{C}$, $B\colon E\longrightarrow\mathbb{C}$, $C\colon\mathbb{C}\longrightarrow E$, and $D\colon E\longrightarrow E$.

The interpolation vectors satisfy 

$$ U\begin{pmatrix}1\\ z_{i}v_{i} \end{pmatrix}=\begin{pmatrix}w_{i}\\ v_{i} \end{pmatrix}. $$

Comparing the two coordinates gives $A+z_{i}Bv_{i}=w_{i}$ and $C+z_{i}Dv_{i}=v_{i}$.

Since $U$ is unitary, $D$ is contractive. Thus $I-zD$ is invertible for every $z\in\mathbb{D}$. The second identity gives $v_{i}=\left(I-z_{i}D\right)^{-1}C$. Substitution into the first identity gives 

$$ w_{i}=A+z_{i}B\left(I-z_{i}D\right)^{-1}C. $$

We are therefore led to define 

$$ f\left(z\right)=A+zB\left(I-zD\right)^{-1}C. $$

This function is holomorphic on $\mathbb{D}$, and the preceding calculation gives $f\left(z_{i}\right)=w_{i}$.

It remains to verify that $f$ belongs to the Schur class. Put $\gamma\left(z\right)=\left(I-zD\right)^{-1}C$. The identity $C+zD\gamma\left(z\right)=\gamma\left(z\right)$ shows that 

$$ U\begin{pmatrix}1\\ z\gamma\left(z\right) \end{pmatrix}=\begin{pmatrix}f\left(z\right)\\ \gamma\left(z\right) \end{pmatrix}. $$

Since $U$ is unitary, 

$$ \begin{aligned}1+\overline{z}w\left\langle \gamma\left(z\right),\gamma\left(w\right)\right\rangle  & =\overline{f\left(z\right)}f\left(w\right)+\left\langle \gamma\left(z\right),\gamma\left(w\right)\right\rangle .\end{aligned} $$

Setting $w=z$ gives $1-\left|f\left(z\right)\right|^{2}=\left(1-\left|z\right|^{2}\right)\left\Vert \gamma\left(z\right)\right\Vert ^{2}$. Hence $\left|f\left(z\right)\right|\leq1$ for every $z\in\mathbb{D}$.

The same calculation gives the kernel factorization 

$$ \frac{1-f\left(z\right)\overline{f\left(w\right)}}{1-z\overline{w}}=\left\langle \gamma\left(w\right),\gamma\left(z\right)\right\rangle . $$

In particular, the de Branges-Rovnyak kernel of the constructed interpolant is positive definite. For any finite collection $u_{1},\ldots,u_{m}\in\mathbb{D}$ and scalars $c_{1},\ldots,c_{m}$, one has 

$$ \begin{aligned}\sum^{m}_{i,j=1}\overline{c_{i}}c_{j}c_{f}\left(u_{i},u_{j}\right) & =\left\Vert \sum^{m}_{i=1}\overline{c_{i}}\gamma\left(u_{i}\right)\right\Vert ^{2}\geq0.\end{aligned} $$

---

## The two-point case

Suppose that $z_{1},z_{2},w_{1},w_{2}\in\mathbb{D}$. When $n=2$, positivity of the Pick matrix reduces to the Schwarz-Pick inequality.

Define the pseudohyperbolic distance by 

$$ \rho\left(z,w\right)=\left|\frac{z-w}{1-\overline{w}z}\right|. $$

For two interpolation conditions $f\left(z_{1}\right)=w_{1}$ and $f\left(z_{2}\right)=w_{2}$, the Pick matrix is 

$$ \begin{pmatrix}\dfrac{1-\left|w_{1}\right|^{2}}{1-\left|z_{1}\right|^{2}} & \dfrac{1-w_{1}\overline{w_{2}}}{1-z_{1}\overline{z_{2}}}\\[4mm] \dfrac{1-w_{2}\overline{w_{1}}}{1-z_{2}\overline{z_{1}}} & \dfrac{1-\left|w_{2}\right|^{2}}{1-\left|z_{2}\right|^{2}} \end{pmatrix}. $$

Its determinant is nonnegative if and only if $\rho\left(w_{1},w_{2}\right)\leq\rho\left(z_{1},z_{2}\right)$.

If one of the prescribed values lies on the unit circle, positivity of the Pick matrix forces both prescribed values to be equal to the same unimodular constant. The corresponding interpolant is constant.

Thus the two-point Nevanlinna-Pick theorem says that holomorphic self-maps of the disk contract the pseudohyperbolic distance. The general finite interpolation theorem may be viewed as the many-point version of the same kernel inequality.

---

## What is special about the Szegő kernel

Let $\mathcal{H}$ be an arbitrary RKHS with kernel $K$. If $f$ is a contractive multiplier of $\mathcal{H}$, then 

$$ \left[\left(1-w_{i}\overline{w_{j}}\right)K\left(z_{i},z_{j}\right)\right]^{n}_{i,j=1}\geq0 $$

whenever $f\left(z_{i}\right)=w_{i}$. This follows from the same calculation on the span of the kernel functions.

For a general kernel, this positivity condition need not be sufficient. The Szegő kernel has the Pick property. Positivity gives a global contractive multiplier.

The corresponding statement for matrix-valued interpolation is the complete Pick property. From this, the classical Nevanlinna-Pick theorem is the basic model for interpolation in complete Pick spaces.

---

## A final kernel factorization

The de Branges-Rovnyak kernel also appears in the factorization 

$$ k\left(z,w\right)=k\left(f\left(z\right),f\left(w\right)\right)c_{f}\left(z,w\right). $$

If $c_{f,w}=c_{f}\left(\cdot,w\right)$, this identity gives an isometry $V_{f}\colon H^{2}\longrightarrow H^{2}\otimes\mathcal{H}\left(f\right)$ defined on the Szegő kernel functions by 

$$ V_{f}k_{w}=k_{f\left(w\right)}\otimes c_{f,w}. $$

Indeed, 

$$ \begin{aligned}\left\langle V_{f}k_{z},V_{f}k_{w}\right\rangle  & =k\left(f\left(z\right),f\left(w\right)\right)c_{f}\left(z,w\right)=k\left(z,w\right).\end{aligned} $$

The Nevanlinna-Pick data prescribe a finite piece of the second kernel in this tensor factorization. The Pick matrix is the proposed Gram matrix of the vectors $c_{f,z_{1}},\ldots,c_{f,z_{n}}$ in the unknown space $\mathcal{H}\left(f\right)$.

The theorem says that positivity of this finite Gram matrix is enough to produce the Schur function, the de Branges-Rovnyak space, and the global multiplier.

This is why the RKHS formulation makes the classical theorem transparent. The interpolation values prescribe an operator on kernel vectors. Contractivity produces the Pick matrix. The same matrix is the finite restriction of the de Branges-Rovnyak kernel. The lurking isometry then turns that finite Hilbert space data into a holomorphic interpolant.
