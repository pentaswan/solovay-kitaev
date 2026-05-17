import { basicApproximation, GCDecompose, transformSU2 } from "./gcdecompose";
import { gateSets } from "./util/availableGates";
import { Matrix } from "./util/matrix";
import { MatrixSet } from "./util/matrixSet";
import './style.css';
import './test.ts';
import { ComplexNumber } from "./util/complex.ts";

export async function SolovayKitaev(U: Matrix, n: number=5): Promise<MatrixSet> {
    let [gate_mat, phase] = await transformSU2(U.clone());
    
    let decomposition = (await basicApproximation(gate_mat))[0];

    for(let i = 0; i < n; i++) {
        decomposition = (await _SolovayKitaev(gate_mat, i + 1, decomposition));
    }

    return decomposition.mulPhase(-phase);
}

async function _SolovayKitaev(U: Matrix, n: number, U_n1_ids: MatrixSet): Promise<MatrixSet> {
    if(n === 0) {
        const decomposition = (await basicApproximation(U))[0];

        return decomposition;
    }
    
    const U_n1_mat = await idsToMat(U_n1_ids);
    let [V, W] = await GCDecompose(U.clone().mul(U_n1_mat.clone().dagger()));

    let V_n1_ids, W_n1_ids;
    for(let i = 0; i < 2; i++) {
        const C_n_ids = [V, W][i];

        let C_n1_ids = new MatrixSet();

        for(let j = 0; j < n; j++) {
            C_n1_ids = (await _SolovayKitaev(C_n_ids.computed, j, C_n1_ids.clone().reverse()));
        }

        if(i === 0) V_n1_ids = C_n1_ids;
        else W_n1_ids = C_n1_ids;
    }

    // compute VWV†W†U
    const newIds = V_n1_ids!.clone().mul(W_n1_ids!.clone()).mul(V_n1_ids!.clone().dagger()).mul(W_n1_ids!.clone().dagger()).mul(U_n1_ids);

    return newIds;
}

// Ids and mat differ by a phase. This function transforms by this phase.
async function idsToMat(ids: MatrixSet): Promise<Matrix> {
    return (await transformSU2(ids.clone().reverse().computed.clone().transpose()))[0].clone();
}

export function groupCommutator(U: MatrixSet, V: MatrixSet): MatrixSet {
    return U.clone().mul(V.clone()).mul(U.clone().dagger()).mul(V.clone().dagger());
}

export function commutator(U: Matrix, V: Matrix): Matrix {
    return U.clone().mul(V.clone()).sub(V.clone().mul(U.clone()));
}