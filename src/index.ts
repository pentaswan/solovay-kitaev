import { basicApproximation, GCDecompose, transformSU2 } from "./gcdecompose";
import { gateSets } from "./util/availableGates";
import { Matrix } from "./util/matrix";
import { MatrixSet } from "./util/matrixSet";
import './style.css';
import './test.ts';
import { ComplexNumber } from "./util/complex.ts";

export async function SolovayKitaev(U: MatrixSet, n: number=5, eps: number=1e-5): Promise<MatrixSet> {
    let [gate_mat, phase] = await transformSU2(U.clone());
    
    let [decomposition, U_prime] = await basicApproximation(gate_mat);

    for(let i = 0; i < n; i++) {
        if(U.computed.distance(U_prime) < eps) {
            console.warn('breaking!');
            break;
        }

        [decomposition, U_prime] = await _SolovayKitaev(gate_mat, i + 1, decomposition, U_prime);
    }

    return decomposition.reverse().mulPhase(-phase);
}

async function _SolovayKitaev(U: MatrixSet, n: number, U_n1_ids: MatrixSet, U_n1_mat: Matrix): Promise<[MatrixSet, Matrix]> {
    if(n === 0) {
        const [decomposition, u_prime] = await basicApproximation(U);

        return [decomposition, u_prime];
    }
    
    let [V, W] = await GCDecompose(new MatrixSet(U, U_n1_mat.clone().dagger()));

    let V_n1_ids, W_n1_ids, V_n1_mat, W_n1_mat;
    for(let i = 0; i < 2; i++) {
        const C_n_ids = [V, W][i];

        let C_n1_ids = new MatrixSet();
        let C_n1_mat = new Matrix();

        for(let j = 0; j < n; j++) {
            [C_n1_ids, C_n1_mat] = await _SolovayKitaev(C_n_ids, j, C_n1_ids, C_n1_mat);
        }

        if(i === 0) {
            V_n1_ids = C_n1_ids;
            V_n1_mat = C_n1_mat;
        } else {
            W_n1_ids = C_n1_ids;
            W_n1_mat = C_n1_mat;
        }
    }

    // compute VWV†W†U and its transpose
    const newIds = U_n1_ids.clone().mul(W_n1_ids!.clone().dagger()).mul(V_n1_ids!.clone().dagger()).mul(W_n1_ids!).mul(V_n1_ids!);
    const newMat = V_n1_mat!.clone().mul(W_n1_mat!).mul(V_n1_mat!.clone().dagger()).mul(W_n1_mat!.clone().dagger()).mul(U_n1_mat);

    return [newIds, newMat];
}

export function groupCommutator(U: MatrixSet, V: MatrixSet): MatrixSet {
    return U.clone().mul(V.clone()).mul(U.clone().dagger()).mul(V.clone().dagger());
}

export function commutator(U: Matrix, V: Matrix): Matrix {
    return U.clone().mul(V.clone()).sub(V.clone().mul(U.clone()));
}