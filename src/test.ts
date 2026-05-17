import { commutator, groupCommutator, SolovayKitaev } from ".";
import { GCDecompose, transformSU2 } from "./gcdecompose";
import { c, gatesAndNames, nameToGate } from "./util/availableGates";
import { ComplexNumber } from "./util/complex";
import { Matrix, NoiseType } from "./util/matrix";
import { MatrixSet } from "./util/matrixSet";


// Ideas:
// - zmq certain parts of the program to verify.
// - perhaps try to simplify the python step by step until it doesn't work / matches the js implementation?
// - try changing distance metric?
// - run python with d = 1, and compare gates generated / trace vs what happens in js. My guess is the error has something to do with different gates being picked (distance metric / spatial tree).

// python, n=0

// [[ 0.90657609-0.17796964j  0.31765609+0.21340389j]
//  [-0.37551611+0.07371744j  0.76688964+0.51520256j]]
// [T(0), H(0), S(0), H(0), T(0), H(0), S(0), GlobalPhase(array(0.58654373), wires=[])]

// js, n=0
// [
//     0.905-0.217i, 0.28+0.235i
//     -0.349+0.108i, 0.751+0.551i
// ]

// Js, according to python server, is trying to basicApprox:
// [[ 0.905-0.217j  0.28 +0.235j]
//  [-0.349+0.108j  0.751+0.551j]]
// Python responded with [S(0), H(0), S(0), T(0), H(0), T(0), H(0)]
// JS Got S, H, S, T, H, T, H




// js, n=1
// calling _s_k, n+1 = 1 (n = 0)

// U
// [
//     0.844-0.392i, 0.321+0.175i
//     -0.321+0.175i, 0.845+0.392i
// ]

// U_approx == decomp product
// [
//     0.854+0.354i, -0.354-0.146i
//     0.146+0.354i, 0.354+0.854i
// ]

// python, n=1
// depth+1 = 1 (n = 0)
// decomposition [T(0), H(0), S(0), H(0), T(0), H(0), S(0)]
// decomp product [[ 0.85355339+0.35355339j -0.35355339-0.14644661j]
//  [ 0.85355339+0.35355339j -0.35355339-0.14644661j]]
// u_prime [[ 0.85355339-0.35355339j  0.35355339+0.14644661j]
//  [-0.35355339+0.14644661j  0.85355339+0.35355339j]]

// const det = U.determinant();

// // Original, working matrix
// const U = new Matrix(
//     c(0.905, -0.217), c(0.28, 0.235),
//     c(-0.349, 0.108), c(0.751, 0.551)
// );

// const U = new Matrix(
//     c(0.439, -0.288), c(0.793, -0.309),
//     c(-0.742, 0.417), c(0.501, -0.156)
// );
/*
const U = new Matrix().addNoise(1, NoiseType.Unitary);

console.log({isUnitary: U.isUnitary(), determinant: U.determinant()});

(async() => {
    const result = await SolovayKitaev(new MatrixSet(U), 3);
    console.log('final result:');
    result.log();
    console.log({computedValid: result.checkComputed()});
    // console.log('U (transformed)');
    // (await transformSU2(new MatrixSet(U))).computed.log();
    console.log('U');
    U.log();
})();
*/


const U = new Matrix().addNoise(1, NoiseType.Unitary);

(async() => {
    const result = await SolovayKitaev(U, 3);

    console.log('Approximating gate:');
    U.log();

    console.log('Result:');
    result.computed.log();
    
    result.log();
})();

// // SolovayKitaev(new MatrixSet(nameToGate.H)).log();
// // // GCDecompose
// const A = new Matrix();
// A.addNoise(0.01, NoiseType.Unitary);
// console.log('A');
// A.log();

// const B = new Matrix();
// B.addNoise(0.01, NoiseType.Unitary);
// console.log('B');
// B.log();

// const delta = new Matrix();
// delta.addNoise(0.01, NoiseType.Unitary);

// const E = new Matrix();
// E.addNoise(0.01, NoiseType.Unitary);

// const [U, V] = GCDecompose(new MatrixSet(E));

// E.sub(nameToGate.I);

// console.log('unitary checks:', {E: E.isUnitary()}, {U: U.computed.isUnitary()}, {V: V.computed.isUnitary()});

// const commutator_uv = commutator(U.computed, V.computed);
// console.log('commutator');
// commutator_uv.log();
// console.log('E');
// E.log();

// console.log('distance', E.distance(commutator_uv));

// const U = A.clone().scaleComplex(c(0, -1)).exp();
// const V = B.clone().scaleComplex(c(0, -1)).exp();

// // Scenario: Delta is our error. Let's assume we've found matricies A and B s.t. delta = e^(-i[A, B])
// const delta = commutator(A, B).scaleComplex(c(0, -1)).exp();

// GCDecompose(new MatrixSet(E))

// e^(-[A, B])

// 
// GCDecompose()

// [U, V]_gp ≈ e^(-[A, B])



// console.log({isUnitary: A.isUnitary()});

// const [V, W] = GCDecompose(new MatrixSet(A));

// groupCommutator(V, W).log();
// // console.log(GCDecompose(new MatrixSet(nameToGate.H)));
// // GCDecompose(new MatrixSet(nameToGate.H)).map(a => a.log());
// // groupCommutator(...GCDecompose(new MatrixSet(nameToGate.H))).log();

// // 1) Eigenvalues
// for(let i = 0; i < gatesAndNames.length; i++) {
//     const [gate, name] = gatesAndNames[i];
//     gate.log();
//     gate.eigenvalues().forEach(a => a.log());
// }

// // 2) Eigenvectors
// for(let i = 0; i < gatesAndNames.length; i++) {
//     const [gate, name] = gatesAndNames[i];

//     console.log('---');
//     console.log({name});
//     gate.log();
//     console.log('eigenvalues');
//     gate.eigenvalues().forEach(a => a.log());
//     console.log('eigenvectors');
//     gate.eigenvectors().log();
//     console.log('---');
// }

// // 3) Group Commutator
// // Idea of the group commutator: If A, B are within a distance eps of the identity, and U = e^(-iA), V = e^(-iB), groupCommutator(U, V) will be O(eps^3) of e^(-[A, B])
// // Let's test this.
// const A = new Matrix();
// A.addNoise(0.01);

// const B = new Matrix();
// B.addNoise(0.01);

// const U = A.clone().scaleComplex(c(0, -1)).exp();
// const V = B.clone().scaleComplex(c(0, -1)).exp();

// console.log('A', {isHermitian: A.isHermitian()});
// A.log();
// console.log('A Distance', A.distance(nameToGate.I));
// console.log('B', {isHermitian: B.isHermitian()});
// B.log();
// console.log('B Distance', B.distance(nameToGate.I));

// // console.log('U');
// // U.log();
// // console.log('V');
// // V.log();

// const g_commutator = groupCommutator(new MatrixSet(U), new MatrixSet(V)).computed;

// const c_commutator = commutator(A, B).scale(-1);

// console.log('C');
// c_commutator.log();

// console.log('GC');
// g_commutator.log();
// console.log('GC Distance', g_commutator.distance(c_commutator.clone().exp()));