"use strict";

describe("TrivialConnections", function() {
	let polygonSoup = MeshIO.readOBJ(solution);
	let mesh = new Mesh();
	mesh.build(polygonSoup);
	let geometry = new Geometry(mesh, polygonSoup["v"], false);
	let E = mesh.edges.length;
	let trivialConnections, singularity, deltaBeta, gamma;

	describe("initialization", function() {
		it("harmonic basis is correct", function() {
                        let myTrivialConnections = new TrivialConnections(geometry);
                        let bases = myTrivialConnections.bases;
                        let hodgeDecomposition = new HodgeDecomposition(geometry);

                        // basis has size 2g
                        let g = 1 - mesh.eulerCharacteristic() / 2;
                        chai.assert.strictEqual(bases.length, 2 * g, 'basis size should be 2g');

                        let baseLen = bases.length;
                        let rankMatrix = DenseMatrix.zeros(baseLen, baseLen);
                        for (let i = 0; i < baseLen; i++) {
                                for (let j = 0; j < baseLen; j++) {
                                        rankMatrix.set(bases[i].transpose().timesDense(bases[j]).get(0, 0), i, j);
                                }
                        }

                        chai.assert.isTrue(rankMatrix.rank() === baseLen, 'basis should be linearly independent');

                        for (let gamma of bases) {
                            let dgamma = hodgeDecomposition.d1.timesDense(gamma);
                            chai.assert.approximately(dgamma.norm(2), 0, 1e-4, 'basis should be closed')

                            let delgamma = hodgeDecomposition.d0T.timesDense(hodgeDecomposition.hodge1.timesDense(gamma));
                            chai.assert.approximately(delgamma.norm(2), 0, 1e-4, 'basis should be exact')
                        }

                        // undo initialization
                        mesh = new Mesh();
                        mesh.build(polygonSoup);
                        geometry = new Geometry(mesh, polygonSoup["v"], false);
		});
	});

	describe("computeCoExactComponent", function() {
		it("computes the dual 0-form potential", function() {
			let loadSingularitiesAndCoExactComponent = function() {
				let singularity = {};
				let deltaBeta = DenseMatrix.zeros(E, 1);

				let v = 0;
				let e = 0;
				let lines = solution.split("\n");
				for (let line of lines) {
					line = line.trim();
					let tokens = line.split(" ");
					let identifier = tokens[0].trim();

					if (identifier === "singularity") {
						singularity[v] = parseFloat(tokens[1]);
						v++;

					} else if (identifier === "deltaBeta") {
						deltaBeta.set(parseFloat(tokens[1]), e, 0);
						e++;
					}
				}

				return [singularity, deltaBeta];
			}

			let [singularity, deltaBeta_sol] = loadSingularitiesAndCoExactComponent();
			trivialConnections = new TrivialConnections(geometry);
			deltaBeta = trivialConnections.computeCoExactComponent(singularity);

			chai.assert.approximately(deltaBeta.minus(deltaBeta_sol).norm(), 0, 1e-6, true);

			let exceptList = [trivialConnections.P, trivialConnections.A,
				trivialConnections.hodge1, trivialConnections.d0, deltaBeta
			];
			exceptList = exceptList.concat(trivialConnections.bases);
			memoryManager.deleteExcept(exceptList);
		});
	});

	describe("buildPeriodMatrix", function() {
		it("computes the dual 0-form potential", function() {
			let loadSingularitiesAndCoExactComponent = function() {
				let singularity = {};
				let deltaBeta = DenseMatrix.zeros(E, 1);

				let v = 0;
				let e = 0;
				let lines = solution.split("\n");
				for (let line of lines) {
					line = line.trim();
					let tokens = line.split(" ");
					let identifier = tokens[0].trim();

					if (identifier === "singularity") {
						singularity[v] = parseFloat(tokens[1]);
						v++;

					} else if (identifier === "deltaBeta") {
						deltaBeta.set(parseFloat(tokens[1]), e, 0);
						e++;
					}
				}

				return [singularity, deltaBeta];
			}

			let [singularity, deltaBeta_sol] = loadSingularitiesAndCoExactComponent();
			trivialConnections = new TrivialConnections(geometry);
			deltaBeta = trivialConnections.computeCoExactComponent(singularity);

			chai.assert.approximately(deltaBeta.minus(deltaBeta_sol).norm(), 0, 1e-6, true);

			let exceptList = [trivialConnections.P, trivialConnections.A,
				trivialConnections.hodge1, trivialConnections.d0, deltaBeta
			];
			exceptList = exceptList.concat(trivialConnections.bases);
			memoryManager.deleteExcept(exceptList);
		});
	});


	describe("computeHarmonicComponent", function() {
		it("computes the harmonic component", function() {
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
						e++;
					}
				}

				return gamma;
			}

			let gamma_sol = loadHarmonicComponent();
			gamma = trivialConnections.computeHarmonicComponent(deltaBeta);

			chai.assert.strictEqual(gamma.minus(gamma_sol).norm() < 1e-6, true);
			let exceptList = [trivialConnections.P, trivialConnections.A,
				trivialConnections.hodge1, trivialConnections.d0
			];
			exceptList = exceptList.concat(trivialConnections.bases);

			memoryManager.deleteExcept(exceptList);
		});
	});

	describe("computeConnections", function() {
		it("computes the dual 1-form connections", function() {
			let loadConnections = function() {
				let phi = DenseMatrix.zeros(E, 1);
				let singularity = {};

				let e = 0;
				let v = 0;
				let lines = solution.split("\n");
				for (let line of lines) {
					line = line.trim();
					let tokens = line.split(" ");
					let identifier = tokens[0].trim();

					if (identifier === "phi") {
						phi.set(parseFloat(tokens[1]), e, 0);
						e++;
					} else if (identifier === "singularity") {
						singularity[v] = parseFloat(tokens[1]);
						v++;
                                        }

				}

				return [phi, singularity];
			}

			let [phi_sol, singularity] = loadConnections();
			let phi = trivialConnections.computeConnections(singularity);//deltaBeta.plus(gamma);

			chai.assert.strictEqual(phi.minus(phi_sol).norm() < 1e-6, true);
			memoryManager.deleteExcept([]);
		});
	});
	describe("computeConnections on a sphere", function() {
                let spherePolygonSoup = MeshIO.readOBJ(sphere);
                let sphereMesh = new Mesh();
                sphereMesh.build(spherePolygonSoup);
                let sphereGeometry, sphereTrivialConnections, sphereSingularity, d0T, u;

		it("connection form is coclosed", function() {
			let loadSingularitiesAndU = function() {
				let singularity = {};
                                let uDict = {}

				let v = 0;
				let e = 0;
				let lines = sphere.split("\n");
				for (let line of lines) {
					line = line.trim();
					let tokens = line.split(" ");
					let identifier = tokens[0].trim();

					if (identifier === "singularity") {
						singularity[v] = parseFloat(tokens[1]);
                                        } else if (identifier == "u") {
						uDict[v] = parseFloat(tokens[1]);
						v++;
					}
				}
                                let u = DenseMatrix.zeros(v, 1);
                                for (let i = 0; i < v; i++) {
                                    u.set(uDict[i], i, 0);
                                }

				return [singularity, u];
			}
                        sphereGeometry = new Geometry(sphereMesh, spherePolygonSoup["v"], false);
                        sphereTrivialConnections = new TrivialConnections(sphereGeometry);
                        [sphereSingularity, u] = loadSingularitiesAndU();

                        let vertexIndex = indexElements(sphereMesh.vertices);
                        let edgeIndex   = indexElements(sphereMesh.edges);
                        let faceIndex   = indexElements(sphereMesh.faces);
                        let hodge1Inv = DEC.buildHodgeStar1Form(sphereGeometry, edgeIndex).invertDiagonal();
                        let d1  = DEC.buildExteriorDerivative1Form(sphereGeometry, faceIndex, edgeIndex);

			let phi = sphereTrivialConnections.computeConnections(sphereSingularity);
                        let codifferential = d1.timesDense(hodge1Inv.timesDense(phi));

                        chai.assert.approximately(codifferential.norm(2), 0, 1e-4, 'connection form is coclosed');
			memoryManager.deleteExcept([d0T, u]);
		});
		it("dδβ = u", function() {
                        sphereGeometry = new Geometry(sphereMesh, spherePolygonSoup["v"], false);
                        sphereTrivialConnections = new TrivialConnections(sphereGeometry);
			let phi = sphereTrivialConnections.computeConnections(sphereSingularity);

                        let vertexIndex = indexElements(sphereMesh.vertices);
                        let edgeIndex   = indexElements(sphereMesh.edges);
                        let d0T = DEC.buildExteriorDerivative0Form(sphereGeometry, edgeIndex, vertexIndex).transpose();
                        let differential = d0T.timesDense(phi);

                        let err = differential.minus(u);

                        chai.assert.approximately(err.norm(2), 0, 1e-4, 'dδβ = u');
		});
	});
});
