"use strict";

class HeatMethod {
	/**
	 * This class implements the {@link http://cs.cmu.edu/~kmcrane/Projects/HeatMethod/ heat method} to compute geodesic distance
	 * on a surface mesh.
	 * @constructor module:Projects.HeatMethod
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {Object} vertexIndex A dictionary mapping each vertex of the input mesh to a unique index.
	 * @property {module:LinearAlgebra.SparseMatrix} A The laplace matrix of the input mesh.
	 * @property {module:LinearAlgebra.SparseMatrix} F The mean curvature flow operator built on the input mesh.
	 */
	constructor(geometry) {
		this.geometry = geometry;
		this.vertexIndex = indexElements(geometry.mesh.vertices);

		// TODO: build laplace and flow matrices
		this.A = SparseMatrix.identity(1, 1); // placeholder
		this.F = SparseMatrix.identity(1, 1); // placeholder
	}

	/**
	 * Computes the vector field X = -∇u / |∇u|.
	 * @private
	 * @method module:Projects.HeatMethod#computeVectorField
	 * @param {module:LinearAlgebra.DenseMatrix} u A dense vector (i.e., u.nCols() == 1) representing the
	 * heat that is allowed to diffuse on the input mesh for a brief period of time.
	 * @returns {Object} A dictionary mapping each face of the input mesh to a {@link module:LinearAlgebra.Vector Vector}.
	 */
	computeVectorField(u) {
		// TODO

		return {}; // placeholder
	}

	/**
	 * Computes the integrated divergence ∇.X.
	 * @private
	 * @method module:Projects.HeatMethod#computeDivergence
	 * @param {Object} X The vector field -∇u / |∇u| represented by a dictionary
	 * mapping each face of the input mesh to a {@link module:LinearAlgebra.Vector Vector}.
	 * @returns {module:LinearAlgebra.DenseMatrix}
	 */
	computeDivergence(X) {
		// TODO

		return DenseMatrix.zeros(1, 1); // placeholder
	}

	/**
	 * Shifts φ such that its minimum value is zero.
	 * @private
	 * @method module:Projects.HeatMethod#subtractMinimumDistance
	 * @param {module:LinearAlgebra.DenseMatrix} phi The (minimum 0) solution to the poisson equation Δφ = ∇.X.
	 */
	subtractMinimumDistance(phi) {
		let min = Infinity;
		for (let i = 0; i < phi.nRows(); i++) {
			min = Math.min(phi.get(i, 0), min);
		}

		for (let i = 0; i < phi.nRows(); i++) {
			phi.set(phi.get(i, 0) - min, i, 0);
		}
	}

	/**
	 * Computes the geodesic distances φ using the heat method.
	 * @method module:Projects.HeatMethod#compute
	 * @param {module:LinearAlgebra.DenseMatrix} delta A dense vector (i.e., delta.nCols() == 1) containing
	 * heat sources, i.e., u0 = δ(x).
	 * @returns {module:LinearAlgebra.DenseMatrix}
	 */
	compute(delta) {
		// TODO
		let phi = DenseMatrix.zeros(delta.nRows(), 1); // placeholder

		// since φ is unique up to an additive constant, it should
		// be shifted such that the smallest distance is zero
		this.subtractMinimumDistance(phi);

		return phi;
	}
}
