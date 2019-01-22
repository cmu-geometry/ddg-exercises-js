"use strict";

class HarmonicBases {
	/**
	 * This class computes the {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf harmonic bases} of a surface mesh.
	 * @constructor module:Projects.HarmonicBases
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 */
	constructor(geometry) {
		this.geometry = geometry;
	}

	/**
	 * Builds a closed, but not exact, primal 1-form ω.
	 * @private
	 * @method module:Projects.HarmonicBases#buildClosedPrimalOneForm
	 * @param {module:Core.Halfedge[]} generator An array of halfedges representing a
	 * {@link https://en.wikipedia.org/wiki/Homology_(mathematics)#Surfaces homology generator}
	 * of the input mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of the input mesh
	 * to a unique index.
	 * @returns {module:LinearAlgebra.DenseMatrix}
	 */
	buildClosedPrimalOneForm(generator, edgeIndex) {
		// TODO

		return DenseMatrix.zeros(1, 1); // placeholder
	}

	/**
	 * Computes the harmonic bases [γ1, γ2 ... γn] of the input mesh.
	 * @method module:Projects.HarmonicBases#compute
	 * @param {module:Projects.HodgeDecomposition} hodgeDecomposition A hodge decomposition object that
	 * can be used to compute the exact component of the closed, but not exact, primal
	 * 1-form ω.
	 * @returns {module:LinearAlgebra.DenseMatrix[]}
	 */
	compute(hodgeDecomposition) {
		// TODO

		return []; // placeholder
	}
}
