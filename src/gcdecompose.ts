import { groupCommutator } from ".";
import { clamp } from "./util";
import { c, gateSets, nameToGate } from "./util/availableGates";
import type { ComplexNumber } from "./util/complex";
import { Matrix } from "./util/matrix";
import { MatrixSet } from "./util/matrixSet";
import { Vector } from "./util/vector";

async function sendServer(data: any, route: string): Promise<string> {
    return new Promise((res) => {
        fetch(`http://127.0.0.1:5000/${route}`, {method: 'POST', headers: {'content-type': 'application/json'}, body: JSON.stringify(data)}).then(async (d) => {
            res(await d.text());
        })
    })
}

export async function GCDecompose(matrix: MatrixSet, tol: number = 1e-5): Promise<[MatrixSet, MatrixSet]> {
    const quaternion = quaternionTransform(matrix.computed);

    const theta = 2 * Math.acos(clamp(quaternion[0], -1, 1));
    const axis = new Vector(quaternion[1], quaternion[2], quaternion[3]);

    if (closeToZero(axis, tol) && Math.abs(theta % Math.PI) < tol) {
        return [new MatrixSet(), new MatrixSet()];
    }

    // compute this weird angle of rotation that inverses the group commutator
    const phi = 2 * Math.asin(Math.sqrt(Math.sqrt(0.5 - 0.5 * Math.cos(theta / 2))));

    // we have lots of options for v and w, so we restrict them to x and y rotations.
    const v = new MatrixSet(RX(phi));
    const w = new MatrixSet(RY(axis.z > 0 ? Math.PI * 2 - phi : phi));

    // diagonalize u, v, w
    const ud = await eigenvectors(matrix);
    const vwd = await eigenvectors(groupCommutator(v, w));

    const s = ud.clone().mul(vwd.clone().dagger());
    const sdg = s.clone().dagger();

    // v' = s * v * s†, w = s * w * s†, transforming the group commutator by s because inner "s"s cancel
    const v_hat = s.clone().mul(v.clone()).mul(sdg);
    const w_hat = s.clone().mul(w.clone()).mul(sdg);

    return [w_hat, v_hat];
}

// // Old working GCDecompose implementation that uses python:
// export async function GCDecompose(matrix: MatrixSet): Promise<[MatrixSet, MatrixSet]> {
//     const data = matrixSetToData(matrix);

//     const r = JSON.parse(await sendServer(data, 'api')) as number[];

//     const v_hat = new MatrixSet(new Matrix(
//         c(r[0], r[1]), c(r[2], r[3]), 
//         c(r[4], r[5]), c(r[6], r[7]), 
//     ));

//     const w_hat = new MatrixSet(new Matrix(
//         c(r[8], r[9]), c(r[10], r[11]), 
//         c(r[12], r[13]), c(r[14], r[15]), 
//     ))

//     return [v_hat, w_hat];
// }

// const [v_hat, w_hat] = await GCDecompose(new MatrixSet(new Matrix(
//     c(1 / Math.sqrt(3)), c(Math.sqrt(2 / 3)),
//     c(Math.sqrt(2 / 3)), c(-1 / Math.sqrt(3))
// )));

// console.log("v_hat");
// v_hat.log();
// console.log("w_hat")
// w_hat.log();

function matrixSetToData(U: MatrixSet): {real: number, imaginary: number}[] {
    const d = (U.computed as any).data;
    let data = [];
    for(let i = 0; i < d.length; i++) {
        const num = d[i] as ComplexNumber;

        data.push({real: num.real(), imaginary: num.imaginary()});
    }

    return data;
}

async function eigenvectors(U: MatrixSet): Promise<MatrixSet> {
    return U.eigenvectors();
}

// // Old working eigenvectors implementation that uses python:
// function eigenvectors(U: MatrixSet): Promise<MatrixSet> {
//     const data = matrixSetToData(new MatrixSet(U.computed));

//     const r = JSON.parse(await sendServer(data, 'eig')) as number[];

//     const matrix = new Matrix(
//         c(r[0], r[1]), c(r[2], r[3]),
//         c(r[4], r[5]), c(r[6], r[7]),
//     );

//     const eig = U.eigenvectors().computed;

//     if(matrix.distance(eig) > 0.01) {
//         console.log('difference!');
//         console.log('original matrix');
//         U.computed.log();
//         console.log('js');
//         eig.log();
//         console.log('py (correct)');
//         matrix.log();
//         console.log({eig: (eig as any).data, matrix: (matrix as any).data});
//     }

//     return new MatrixSet(matrix);
// }

function _basicApproximation(U: MatrixSet): MatrixSet {
    const target = transformSU2NoPhase(U).computed;

    let closestSet = new MatrixSet();
    let closestDistance = Infinity;
    for(let i = 0; i < gateSets.length; i++) {
        const gateSet = gateSets[i];

        const approx = transformSU2NoPhase(U);

        const distance = approx.computed.distance(target);

        if(distance < closestDistance) {
            closestDistance = distance;
            closestSet = gateSet;
        }
    }

    return closestSet.clone();
}

export async function basicApproximation(U: MatrixSet): Promise<[MatrixSet, number]> {
    const data = matrixSetToData(U);

    const r = JSON.parse(await sendServer(data, 'approx')) as number[];

    const global_phase = r.pop()!;

    let decomposition = new MatrixSet();
    let u_prime = new Matrix();

    for(let i = 0; i < r.length; i += 8) {
        const matrix = new Matrix(
            c(r[i+0], r[i+1]), c(r[i+2], r[i+3]),
            c(r[i+4], r[i+5]), c(r[i+6], r[i+7]),
        );

        if(i === 0) {
            u_prime = matrix;
        } else {
            decomposition.mul(new MatrixSet(matrix));
        }
    }

    // const jsDecomposition = _basicApproximation(U.clone());

    // console.log('correct');
    // transformSU2NoPhase(new MatrixSet(decomposition.computed)).computed.log();
    // console.log('js');
    // transformSU2NoPhase(new MatrixSet(jsDecomposition.computed)).computed.log();

    // console.log(decomposition.computed.clone().transformSU2()[0].distance(u_prime.clone().transformSU2()[0]));

    const [decomp_SU2, phase] = decomposition.computed.clone().transformSU2();
    decomp_SU2.transpose();

    // if(Math.abs(phase - global_phase) > 0.01) {
    //     console.error('---');
    //     console.log({phase, global_phase});
    //     throw new Error('Phase difference!');
    // }

    // if(decomp_SU2.distance(u_prime) > 0.01) {
    //     console.error('---');
    //     console.log('decomposition');
    //     decomp_SU2.log();
    //     console.log('u_prime');
    //     u_prime.log();
    //     throw new Error('decomposition difference!');
    // }

    return [decomposition.clone().reverse(), global_phase/*phase*/];
}

export async function transformSU2(U: MatrixSet): Promise<[MatrixSet, number]> {
    const [_, globalPhase] = await basicApproximation(U);
    const [decomposed, gatePhase] = U.computed.clone().transformSU2();

    const phase = globalPhase - gatePhase;

    return [new MatrixSet(decomposed), phase];
}

// TEMP
function transformSU2NoPhase(U: MatrixSet): MatrixSet {
    const [decomposed, _] = U.computed.clone().transformSU2();

    return new MatrixSet(decomposed);
}

// // Old working SU2 implementation that uses python:
// export async function _transformSU2(U: MatrixSet): Promise<[MatrixSet, number]> {
//     const data = matrixSetToData(U);

//     const r = JSON.parse(await sendServer(data, 'su2')) as number[];

//     const decomposed = new MatrixSet(new Matrix(
//         c(r[0], r[1]), c(r[2], r[3]), 
//         c(r[4], r[5]), c(r[6], r[7]),
//     ))

//     const phase = r[8];

//     return [decomposed, phase];
// }

function closeToZero(vector: Vector, tol: number) {
    for (let i = 0; i < vector.data.length; i++) {
        if (Math.abs(vector.data[i]) > tol) return false;
    }

    return true;
}

// // Old _SK function from index. has both newIds and newMat.
// async function _SolovayKitaev(U: MatrixSet, n: number, U_n1_ids: MatrixSet, U_n1_mat: Matrix): Promise<[MatrixSet, Matrix]> {
//     if(n === 0) {
//         const [decomposition, u_prime] = await basicApproximation(U);

//         return [decomposition, u_prime];
//     }
    
//     let [V, W] = await GCDecompose(new MatrixSet(U, U_n1_mat.clone().dagger()));

//     let V_n1_ids, W_n1_ids, V_n1_mat, W_n1_mat;
//     for(let i = 0; i < 2; i++) {
//         const C_n_ids = [V, W][i];

//         let C_n1_ids = new MatrixSet();
//         let C_n1_mat = new Matrix();

//         for(let j = 0; j < n; j++) {
//             [C_n1_ids, C_n1_mat] = await _SolovayKitaev(C_n_ids, j, C_n1_ids, C_n1_mat);
//         }

//         if(i === 0) {
//             V_n1_ids = C_n1_ids;
//             V_n1_mat = C_n1_mat;
//         } else {
//             W_n1_ids = C_n1_ids;
//             W_n1_mat = C_n1_mat;
//         }
//     }

//     // compute VWV†W†U and its transpose
//     const newIds = U_n1_ids.clone().mul(W_n1_ids!.clone().dagger()).mul(V_n1_ids!.clone().dagger()).mul(W_n1_ids!).mul(V_n1_ids!);
//     const newMat = V_n1_mat!.clone().mul(W_n1_mat!).mul(V_n1_mat!.clone().dagger()).mul(W_n1_mat!.clone().dagger()).mul(U_n1_mat);

//     return [newIds, newMat];
// }

function quaternionTransform(matrix: Matrix): number[] {
    return [
        matrix.getElement(0, 0).real(),
        -matrix.getElement(0, 1).imaginary(),
        -matrix.getElement(0, 1).real(),
        -matrix.getElement(0, 0).imaginary()
    ]
}

function RX(phi: number): Matrix {
    const cos = Math.cos;
    const sin = Math.sin;

    return new Matrix(
        c(cos(phi / 2)), c(0, -sin(phi / 2)),
        c(0, -sin(phi / 2)), c(cos(phi / 2))
    )
}

function RY(phi: number): Matrix {
    const cos = Math.cos;
    const sin = Math.sin;

    return new Matrix(
        c(cos(phi / 2)), c(-sin(phi / 2)),
        c(sin(phi / 2)), c(cos(phi / 2))
    )
}