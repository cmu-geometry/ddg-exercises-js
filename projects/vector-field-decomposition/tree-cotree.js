"use strict";

class TreeCotree {
	/**
	 * This class computes the {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf tree cotree} decomposition of a surface mesh
	 * to build its {@link https://en.wikipedia.org/wiki/Homology_(mathematics)#Surfaces homology generators}.
	 * @constructor module:Projects.TreeCotree
	 * @param {module:Core.Mesh} mesh The input mesh this class acts on.
	 * @property {module:Core.Mesh} mesh The input mesh this class acts on.
	 * @property {vertexParent} vertexParent A dictionary mapping each vertex of the input mesh to
	 * its parent in the primal spanning tree.
	 * @property {faceParent} faceParent A dictionary mapping each face of the input mesh to
	 * its parent in the dual spanning tree.
	 */
	constructor(mesh) {
		this.mesh = mesh;
		this.vertexParent = {};
		this.faceParent = {};
	}

	/**
	 * Builds a primal spanning tree on a boundaryless mesh.
	 * @private
	 * @method module:Projects.TreeCotree#buildPrimalSpanningTree
	 */
	buildPrimalSpanningTree() {
		// TODO
	}

	/**
	 * Checks whether a halfedge is in the primal spanning tree.
	 * @private
	 * @method module:Projects.TreeCotree#inPrimalSpanningTree
	 * @param {module:Core.Halfedge} h A halfedge on the input mesh.
	 * @returns {boolean}
	 */
	inPrimalSpanningTree(h) {
		// TODO

		return false; // placeholder
	}

	/**
	 * Builds a dual spanning tree on a boundaryless mesh.
	 * @private
	 * @method module:Projects.TreeCotree#buildDualSpanningCotree
	 */
	buildDualSpanningCotree() {
		// TODO
	}

	/**
	 * Checks whether a halfedge is in the dual spanning tree.
	 * @private
	 * @method module:Projects.TreeCotree#inDualSpanningTree
	 * @param {module:Core.Halfedge} h A halfedge on the input mesh.
	 * @returns {boolean}
	 */
	inDualSpanningTree(h) {
		// TODO

		return false; // placeholder
	}

	/**
	 * Returns a halfedge lying on the shared edge between face f and g.
	 * @private
	 * @method module:Projects.TreeCotree#sharedHalfedge
	 * @param {module:Core.Face} f A face on the input mesh.
	 * @param {module:Core.Face} g A neighboring face to f on the input mesh.
	 * @returns {module:Core.Halfedge}
	 */
	sharedHalfedge(f, g) {
		for (let h of f.adjacentHalfedges()) {
			if (h.twin.face === g) {
				return h;
			}
		}

		alert("Code should not reach here!");
		return new Halfedge();
	}

	/**
	 * Computes the {@link https://en.wikipedia.org/wiki/Homology_(mathematics)#Surfaces homology generators} of the input mesh and stores them
	 * in the {@link module:Core.Mesh Mesh}'s generators property.
	 * @method module:Projects.TreeCotree#buildGenerators
	 */
	buildGenerators() {
		// build spanning trees
		this.buildPrimalSpanningTree();
		this.buildDualSpanningCotree();

		// TODO
	}
}
