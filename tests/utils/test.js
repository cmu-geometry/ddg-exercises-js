"use strict";

describe("Solvers", function() {
	describe("Residual", function() {
		it("computes the residual Ax - λx", function() {
                        let T = new ComplexTriplet(3, 3);
                        T.addEntry(new Complex(25,  0), 0, 0);
                        T.addEntry(new Complex(-8,  3), 0, 1);
                        T.addEntry(new Complex(-5, -3), 0, 2);
                        T.addEntry(new Complex(-8, -3), 1, 0);
                        T.addEntry(new Complex(25,  0), 1, 1);
                        T.addEntry(new Complex(-5,  3), 1, 2);
                        T.addEntry(new Complex(-5,  3), 2, 0);
                        T.addEntry(new Complex(-5, -3), 2, 1);
                        T.addEntry(new Complex(22,  0), 2, 2);
                        let M = ComplexSparseMatrix.fromTriplet(T);

                        let x = ComplexDenseMatrix.zeros(3, 1);
                        x.set(new Complex(1, 0), 0, 0);
                        x.set(new Complex(0, 1), 1, 0);
                        x.set(new Complex(0, 1), 2, 0);

                        let expected_residual = 10.43498389;

                        let computed_residual = Solvers.residual(M, x);

                        chai.assert.closeTo(expected_residual, computed_residual, 1e-8);
                });
	});

	describe("Inverse Power Method", function() {
		it("computes the least eigenvalue of a matrix M", function() {
                        let T = new ComplexTriplet(3, 3);
                        T.addEntry(new Complex(25,  0), 0, 0);
                        T.addEntry(new Complex(-8,  3), 0, 1);
                        T.addEntry(new Complex(-5, -3), 0, 2);
                        T.addEntry(new Complex(-8, -3), 1, 0);
                        T.addEntry(new Complex(25,  0), 1, 1);
                        T.addEntry(new Complex(-5,  3), 1, 2);
                        T.addEntry(new Complex(-5,  3), 2, 0);
                        T.addEntry(new Complex(-5, -3), 2, 1);
                        T.addEntry(new Complex(22,  0), 2, 2);
                        let M = ComplexSparseMatrix.fromTriplet(T);

                        let computed_evec = Solvers.solveInversePowerMethod(M);

                        let expected_evec = ComplexDenseMatrix.zeros(3, 1);
                        expected_evec.set(new Complex(-1,  1), 0, 0);
                        expected_evec.set(new Complex(-1, -1), 1, 0);
                        expected_evec.set(new Complex( 2,  0), 2, 0);
                        expected_evec.scaleBy(new Complex(1 / expected_evec.norm(2), 0));

                        let dot = expected_evec.conjugate().transpose().timesDense(computed_evec).get(0,0);

                        chai.assert.closeTo(1, computed_evec.norm(2), 1e-8);
                        chai.assert.closeTo(1, dot.norm(2), 1e-8);
                });
	});
});
