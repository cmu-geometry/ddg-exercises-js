"use strict";

describe("VectorFieldDecomposition", function() {
	let polygonSoup = MeshIO.readOBJ(solution);
	let mesh = new Mesh();
	mesh.build(polygonSoup);
	let E = mesh.edges.length;
	let geometry = new Geometry(mesh, polygonSoup["v"], false);
	let hodgeDecomposition, omega, dAlpha, deltaBeta;

	describe("HodgeDecomposition: initialization", function() {
                let tet_soup = MeshIO.readOBJ(`v 0 0 0
                                               v 1 0 0
                                               v 0 1 0
                                               v 0 0 1
                                               f 1 2 3
                                               f 1 4 2
                                               f 2 4 3
                                               f 1 3 4`)
                let tet_mesh = new Mesh();
                tet_mesh.build(tet_soup);
                let tet_geometry = new Geometry(tet_mesh, tet_soup["v"], false);

		it("A matrix", function() {
                    let tet_hodge_decomposition = new HodgeDecomposition(tet_geometry);

                    let true_A_triples = new Triplet(4, 4);
                    true_A_triples.addEntry(3.00000001, 0, 0);
                    true_A_triples.addEntry(-1, 1, 0);
                    true_A_triples.addEntry(-1, 0, 1);
                    true_A_triples.addEntry(-1, 2, 0);
                    true_A_triples.addEntry(-1, 0, 2);
                    true_A_triples.addEntry(-1, 3, 0);
                    true_A_triples.addEntry(-1, 0, 3);

                    true_A_triples.addEntry( 1.5773502791896257, 1, 1);
                    true_A_triples.addEntry(-0.2886751345948128, 1, 2);
                    true_A_triples.addEntry(-0.2886751345948128, 2, 1);
                    true_A_triples.addEntry(-0.2886751345948128, 1, 3);
                    true_A_triples.addEntry(-0.2886751345948128, 3, 1);

                    true_A_triples.addEntry( 1.5773502791896257, 2, 2);
                    true_A_triples.addEntry(-0.2886751345948128, 2, 3);
                    true_A_triples.addEntry(-0.2886751345948128, 3, 2);

                    true_A_triples.addEntry( 1.5773502791896257, 3, 3);

                    let true_A = SparseMatrix.fromTriplet(true_A_triples);

                    let err = true_A.minus(tet_hodge_decomposition.A)
                    chai.assert.approximately(err.frobeniusNorm(2), 0, 1e-4);
		});
		it("B matrix", function() {
                    let tet_hodge_decomposition = new HodgeDecomposition(tet_geometry);

                    let true_B_triples = new Triplet(4, 4);
                    true_B_triples.addEntry(5.464101615137755, 0, 0);
                    true_B_triples.addEntry(-1, 1, 0);
                    true_B_triples.addEntry(-1, 0, 1);
                    true_B_triples.addEntry(-3.464101615137755, 2, 0);
                    true_B_triples.addEntry(-3.464101615137755, 0, 2);
                    true_B_triples.addEntry(-1, 3, 0);
                    true_B_triples.addEntry(-1, 0, 3);

                    true_B_triples.addEntry(5.464101615137755,  1, 1);
                    true_B_triples.addEntry(-3.464101615137755, 1, 2);
                    true_B_triples.addEntry(-3.464101615137755, 2, 1);
                    true_B_triples.addEntry(-1, 1, 3);
                    true_B_triples.addEntry(-1, 3, 1);

                    true_B_triples.addEntry(10.392304845413264, 2, 2);
                    true_B_triples.addEntry(-3.464101615137755, 2, 3);
                    true_B_triples.addEntry(-3.464101615137755, 3, 2);

                    true_B_triples.addEntry(5.464101615137755,  3, 3);

                    let true_B = SparseMatrix.fromTriplet(true_B_triples);

                    let err = true_B.minus(tet_hodge_decomposition.B)
                    chai.assert.approximately(err.frobeniusNorm(2), 0, 1e-4);
		});
	});

	describe("HodgeDecomposition: computeExactComponent", function() {
		it("exact component is closed", function() {
			let my_omega = DenseMatrix.random(geometry.mesh.edges.length, 1);
			hodgeDecomposition = new HodgeDecomposition(geometry);
			let my_exact_component = hodgeDecomposition.computeExactComponent(my_omega);

                        let d = hodgeDecomposition.d1.timesDense(my_exact_component);

                        chai.assert.approximately(d.norm(2), 0, 1e-4, 'dAlpha approximately closed')
			memoryManager.deleteExcept([ hodgeDecomposition.d1,
				hodgeDecomposition.hodge1Inv, hodgeDecomposition.d1T, hodgeDecomposition.B,
                                hodgeDecomposition.d0T, hodgeDecomposition.hodge1
			]);
		});
		it("computes the exact component of a 1-form", function() {
			let loadOmegaAndExactComponent = function() {
				let omega = DenseMatrix.zeros(E, 1);
				let dAlpha = DenseMatrix.zeros(E, 1);

				let e = 0;
				let lines = solution.split("\n");
				for (let line of lines) {
					line = line.trim();
					let tokens = line.split(" ");
					let identifier = tokens[0].trim();

					if (identifier === "omega") {
						omega.set(parseFloat(tokens[1]), e, 0);

					} else if (identifier === "dAlpha") {
						dAlpha.set(parseFloat(tokens[1]), e, 0);
						e++;
					}
				}

				return [omega, dAlpha];
			}

			let oneForms = loadOmegaAndExactComponent();
			omega = oneForms[0];
			let dAlpha_sol = oneForms[1];
			hodgeDecomposition = new HodgeDecomposition(geometry);
			dAlpha = hodgeDecomposition.computeExactComponent(omega);

			chai.assert.strictEqual(dAlpha.minus(dAlpha_sol).norm() < 1e-6, true);
			memoryManager.deleteExcept([omega, dAlpha, hodgeDecomposition.d1,
				hodgeDecomposition.hodge1Inv, hodgeDecomposition.d1T, hodgeDecomposition.B,
                                hodgeDecomposition.d0T, hodgeDecomposition.hodge1
			]);
		});
	});

	describe("HodgeDecomposition: computeCoExactComponent", function() {
		it("coexact component is coclosed", function() {
			let my_omega = DenseMatrix.random(geometry.mesh.edges.length, 1);
			hodgeDecomposition = new HodgeDecomposition(geometry);
			let my_coexact_component = hodgeDecomposition.computeCoExactComponent(my_omega);

                        let d = hodgeDecomposition.d0T.timesDense(hodgeDecomposition.hodge1.timesDense(my_coexact_component));

                        chai.assert.approximately(d.norm(2), 0, 1e-4, 'deltaBeta approximately coclosed')
			memoryManager.deleteExcept([omega, dAlpha, hodgeDecomposition.d1,
				hodgeDecomposition.hodge1Inv, hodgeDecomposition.d1T, hodgeDecomposition.B,
                                hodgeDecomposition.d0T, hodgeDecomposition.hodge1
			]);
		});
		it("computes the coexact component of a 1-form", function() {
			let loadCoExactComponent = function() {
				let deltaBeta = DenseMatrix.zeros(E, 1);

				let e = 0;
				let lines = solution.split("\n");
				for (let line of lines) {
					line = line.trim();
					let tokens = line.split(" ");
					let identifier = tokens[0].trim();

					if (identifier === "deltaBeta") {
						deltaBeta.set(parseFloat(tokens[1]), e, 0);
						e++
					}
				}

				return deltaBeta;
			}

			let deltaBeta_sol = loadCoExactComponent();
			deltaBeta = hodgeDecomposition.computeCoExactComponent(omega);

			chai.assert.strictEqual(deltaBeta.minus(deltaBeta_sol).norm() < 1e-6, true);
			memoryManager.deleteExcept([omega, dAlpha, deltaBeta, hodgeDecomposition.d1,
				hodgeDecomposition.hodge1Inv, hodgeDecomposition.d1T, hodgeDecomposition.B,
                                hodgeDecomposition.d0T, hodgeDecomposition.hodge1
			]);
		});
	});

	describe("HodgeDecomposition: computeHarmonicComponent", function() {
		it("harmonic component is closed", function() {
			let my_omega = DenseMatrix.random(geometry.mesh.edges.length, 1);
			hodgeDecomposition = new HodgeDecomposition(geometry);
			let my_exact_component   = hodgeDecomposition.computeExactComponent(my_omega);
			let my_coexact_component = hodgeDecomposition.computeCoExactComponent(my_omega);
			let gamma = hodgeDecomposition.computeHarmonicComponent(my_omega, my_exact_component, my_coexact_component);

                        let d = hodgeDecomposition.d1.timesDense(gamma);

                        chai.assert.approximately(d.norm(2), 0, 1e-4, 'gamma approximately closed')
			memoryManager.deleteExcept([omega, dAlpha, deltaBeta, hodgeDecomposition.d1,
				hodgeDecomposition.hodge1Inv, hodgeDecomposition.d1T, hodgeDecomposition.B,
                                hodgeDecomposition.d0T, hodgeDecomposition.hodge1
			]);
		});
		it("harmonic component is coclosed", function() {
			let my_omega = DenseMatrix.random(geometry.mesh.edges.length, 1);
			hodgeDecomposition = new HodgeDecomposition(geometry);
			let my_exact_component   = hodgeDecomposition.computeExactComponent(my_omega);
			let my_coexact_component = hodgeDecomposition.computeCoExactComponent(my_omega);
			let gamma = hodgeDecomposition.computeHarmonicComponent(my_omega, my_exact_component, my_coexact_component);

                        let d = hodgeDecomposition.d0T.timesDense(hodgeDecomposition.hodge1.timesDense(gamma));

                        chai.assert.approximately(d.norm(2), 0, 1e-4, 'gamma approximately coclosed')
			memoryManager.deleteExcept([omega, dAlpha, deltaBeta, hodgeDecomposition.d1,
				hodgeDecomposition.hodge1Inv, hodgeDecomposition.d1T, hodgeDecomposition.B,
                                hodgeDecomposition.d0T, hodgeDecomposition.hodge1
			]);
		});
		it("decomposition sums to original 1-form", function() {
			let my_omega = DenseMatrix.random(geometry.mesh.edges.length, 1);
			hodgeDecomposition = new HodgeDecomposition(geometry);
			let my_exact_component   = hodgeDecomposition.computeExactComponent(my_omega);
			let my_coexact_component = hodgeDecomposition.computeCoExactComponent(my_omega);
			let gamma = hodgeDecomposition.computeHarmonicComponent(my_omega, my_exact_component, my_coexact_component);

                        let reconstructed_omega = gamma.plus(my_exact_component.plus(my_coexact_component));

                        let err = my_omega.minus(reconstructed_omega);

                        chai.assert.approximately(err.norm(2), 0, 1e-4, 'decomposition sums to omega')
			memoryManager.deleteExcept([omega, dAlpha, deltaBeta, hodgeDecomposition.d1,
				hodgeDecomposition.hodge1Inv, hodgeDecomposition.d1T, hodgeDecomposition.B,
                                hodgeDecomposition.d0T, hodgeDecomposition.hodge1
			]);
		});
		it("computes the harmonic component of a 1-form", function() {
			let loadHarmonicComponent = function() {
				let gamma = DenseMatrix.zeros(E, 1);

				let e = 0;
				let lines = solution.split("\n");
				for (let line of lines) {
					line = line.trim();
					let tokens = line.split(" ");
					let identifier = tokens[0].trim();

					if (identifier === "gamma") {
						gamma.set(parseFloat(tokens[1]), e, 0);
						e++
					}
				}

				return gamma;
			}

			let gamma_sol = loadHarmonicComponent();
			let gamma = hodgeDecomposition.computeHarmonicComponent(omega, dAlpha, deltaBeta);

			chai.assert.strictEqual(gamma.minus(gamma_sol).norm() < 1e-6, true);
			memoryManager.deleteExcept([]);
		});
	});

	describe("TreeCotree: buildPrimalSpanningTree", function() {
		it("has exactly one self-edge (the root)", function() {
			let treeCotree = new TreeCotree(mesh);
                        treeCotree.buildPrimalSpanningTree();
                        let n_self_edges = 0;
                        for (let v of mesh.vertices) {
                            if (treeCotree.vertexParent[v] == v) {
                                n_self_edges++;
                            }
                        }

			chai.assert.strictEqual(n_self_edges, 1);
		});
		it("graph is connected", function() {
			let treeCotree = new TreeCotree(mesh);
                        treeCotree.buildPrimalSpanningTree();
                        let root = null;
                        for (let v of mesh.vertices) {
                            if (treeCotree.vertexParent[v] == v) {
                                root = v;
                                break;
                            }
                        }

                        let n = mesh.vertices.length;
                        let connected = true;
                        for (let v of mesh.vertices) {
                            let path_length = 0;
                            while (v != root && path_length < n+1) {
                                v = treeCotree.vertexParent[v];
                                path_length += 1;
                            }
                            if (v != root) {
                                connected = false;
                                break;
                            }
                        }
                        chai.assert.isTrue(connected);
		});
	});

	describe("TreeCotree: buildDualSpanningTree", function() {
		it("has exactly one self-edge (the root)", function() {
			let treeCotree = new TreeCotree(mesh);
                        treeCotree.buildDualSpanningCotree();
                        let n_self_edges = 0;
                        for (let f of mesh.faces) {
                            if (treeCotree.faceParent[f] == f) {
                                n_self_edges++;
                            }
                        }

			chai.assert.strictEqual(n_self_edges, 1);
		});
		it("graph is connected", function() {
			let treeCotree = new TreeCotree(mesh);
                        treeCotree.buildDualSpanningCotree();
                        let root = null;
                        for (let f of mesh.faces) {
                            if (treeCotree.faceParent[f] == f) {
                                root = f;
                                break;
                            }
                        }

                        let n = mesh.faces.length;
                        let connected = true;
                        for (let f of mesh.faces) {
                            let path_length = 0;
                            while (f != root && path_length < n+1) {
                                f = treeCotree.faceParent[f];
                                path_length += 1;
                            }
                            if (f != root) {
                                connected = false;
                                break;
                            }
                        }
                        chai.assert.isTrue(connected);
		});
		it("disjoint from primal spanning tree (uses inPrimalSpanningTree and inDualSpanningTree)", function() {
			let treeCotree = new TreeCotree(mesh);
                        treeCotree.buildPrimalSpanningTree();
                        treeCotree.buildDualSpanningCotree();
                        for (let e of mesh.edges) {
                            let h = e.halfedge;
                            chai.assert.isFalse(treeCotree.inPrimalSpanningTree(h) && treeCotree.inDualSpanningTree(h));
                        }
		});
	});

	describe("TreeCotree: inPrimalSpanningTree", function() {
		it("recognizes positively oriented halfedges in tree", function() {
			let treeCotree = new TreeCotree(mesh);
                        let he = mesh.edges[0].halfedge;
                        let v1 = he.vertex;
                        let v2 = he.twin.vertex;
                        treeCotree.vertexParent[v1] = v2;

			chai.assert.isTrue(treeCotree.inPrimalSpanningTree(he));
		});
		it("recognizes negatively oriented halfedges in tree", function() {
			let treeCotree = new TreeCotree(mesh);
                        let he = mesh.edges[0].halfedge;
                        let v1 = he.twin.vertex;
                        let v2 = he.vertex;
                        treeCotree.vertexParent[v1] = v2;

			chai.assert.isTrue(treeCotree.inPrimalSpanningTree(he));
		});
		it("rejects halfedges not in tree", function() {
			let treeCotree = new TreeCotree(mesh);
                        let he = mesh.edges[0].halfedge;
                        let v1 = he.twin.vertex;
                        let v2 = he.vertex;

			chai.assert.isFalse(treeCotree.inPrimalSpanningTree(he));
		});
	});


	describe("TreeCotree: inDualSpanningTree", function() {
		it("recognizes positively oriented halfedges in tree", function() {
			let treeCotree = new TreeCotree(mesh);
                        let he = mesh.edges[0].halfedge;
                        let f1 = he.face;
                        let f2 = he.twin.face;
                        treeCotree.faceParent[f1] = f2;

			chai.assert.isTrue(treeCotree.inDualSpanningTree(he));
		});
		it("recognizes negatifely oriented halfedges in tree", function() {
			let treeCotree = new TreeCotree(mesh);
                        let he = mesh.edges[0].halfedge;
                        let f1 = he.twin.face;
                        let f2 = he.face;
                        treeCotree.faceParent[f1] = f2;

			chai.assert.isTrue(treeCotree.inDualSpanningTree(he));
		});
		it("rejects halfedges not in tree", function() {
			let treeCotree = new TreeCotree(mesh);
                        let he = mesh.edges[0].halfedge;
                        let f1 = he.twin.face;
                        let f2 = he.face;

			chai.assert.isFalse(treeCotree.inDualSpanningTree(he));
		});
	});

	describe("TreeCotree: buildGenerators", function() {
		it("Has 2g generators", function() {
                        let g = 1 - mesh.eulerCharacteristic() / 2;
			let treeCotree = new TreeCotree(mesh);
                        treeCotree.buildGenerators();
			chai.assert.strictEqual(mesh.generators.length, 2 * g);
		});
	});

	describe("HarmonicBases: buildClosedPrimalOneForm", function() {
		it("Correct nonzero entries", function() {
                        let edgeIndex = indexElements(mesh.edges);
                        let generator = [];
                        for (let i = 0; i < 5; i++) {
                            generator.push(mesh.edges[i].halfedge);
                        }

			let harmonicBases = new HarmonicBases(geometry);
                        let w = harmonicBases.buildClosedPrimalOneForm(generator, edgeIndex);

                        for (let e of mesh.edges) {
                            let i = edgeIndex[e];
                            if (i < 5) {
                                chai.assert.notStrictEqual(w.get(i, 0), 0);
                            } else {
                                chai.assert.strictEqual(w.get(i, 0), 0);
                            }
                        }
		});
		it("Correct orientation", function() {
                        let edgeIndex = indexElements(mesh.edges);
                        let generator = [];
                        for (let i = 0; i < 5; i++) {
                            if (i % 2 == 0) {
                                generator.push(mesh.edges[i].halfedge);
                            } else {
                                generator.push(mesh.edges[i].halfedge.twin);
                            }
                        }

			let harmonicBases = new HarmonicBases(geometry);
                        let w = harmonicBases.buildClosedPrimalOneForm(generator, edgeIndex);

                        for (let e of mesh.edges) {
                            let i = edgeIndex[e];
                            if (i < 5) {
                                if (i % 2 == 0) {
                                    chai.assert.strictEqual(w.get(i, 0), 1);
                                } else {
                                    chai.assert.strictEqual(w.get(i, 0), -1);
                                }
                            } else {
                                chai.assert.strictEqual(w.get(i, 0), 0);
                            }
                        }
		});
	});

	describe("HarmonicBases: compute", function() {
		it("the basis vectors are linearly independent", function() {
			// build generators
			hodgeDecomposition = new HodgeDecomposition(geometry);
			let treeCotree = new TreeCotree(mesh);
			treeCotree.buildGenerators();

			// build harmonic bases
			let harmonicBases = new HarmonicBases(geometry);
			let bases = harmonicBases.compute(hodgeDecomposition);

			let N = bases.length;
			let rankMatrix = DenseMatrix.zeros(N, N);
			for (let i = 0; i < N; i++) {
				for (let j = 0; j < N; j++) {
					rankMatrix.set(bases[i].transpose().timesDense(bases[j]).get(0, 0), i, j);
				}
			}

			chai.assert.strictEqual(rankMatrix.rank() === N, true);
			memoryManager.deleteExcept([]);
		});

		it("the basis vectors are closed", function() {
			// build generators
			hodgeDecomposition = new HodgeDecomposition(geometry);
			let treeCotree = new TreeCotree(mesh);
			treeCotree.buildGenerators();

			// build harmonic bases
			let harmonicBases = new HarmonicBases(geometry);
			let bases = harmonicBases.compute(hodgeDecomposition);

                        for (let gamma of bases) {
                            let dgamma = hodgeDecomposition.d1.timesDense(gamma);
                            chai.assert.approximately(dgamma.norm(2), 0, 1e-4, 'gamma approximately closed')
                        }

			memoryManager.deleteExcept([]);
		});

		it("the basis vectors are coclosed", function() {
			// build generators
			hodgeDecomposition = new HodgeDecomposition(geometry);
			let treeCotree = new TreeCotree(mesh);
			treeCotree.buildGenerators();

			// build harmonic bases
			let harmonicBases = new HarmonicBases(geometry);
			let bases = harmonicBases.compute(hodgeDecomposition);

                        for (let gamma of bases) {
                            let dgamma = hodgeDecomposition.d0T.timesDense(hodgeDecomposition.hodge1.timesDense(gamma));
                            chai.assert.approximately(dgamma.norm(2), 0, 1e-4, 'gamma approximately coclosed')
                        }

			memoryManager.deleteExcept([]);
		});
	});
});
