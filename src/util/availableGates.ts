import { ComplexNumber } from "./complex";
import { Matrix } from "./matrix";
import { MatrixSet } from "./matrixSet";

export const c = (a?: number, b?: number) => { return new ComplexNumber(a, b); }
const sqrt = Math.sqrt;

export const availableGates: Matrix[] = [];
export const gatesAndNames: [Matrix, string][] = [];
export const nameToGate: Record<string, Matrix> = {};
export const gateSets: MatrixSet[] = [];

const I = new Matrix(
    c(1), c(0),
    c(0), c(1)
);

nameToGate['I'] = nameToGate['I_dagger'] = I;
gatesAndNames.push(
    [I, 'I'],
    [I, 'I_dagger']
)

addGate('X', new Matrix(
    c(0), c(1),
    c(1), c(0)
));

addGate('Y', new Matrix(
    c(0), c(0, -1),
    c(0, 1), c(0)
));

addGate('Z', new Matrix(
    c(1), c(0),
    c(0), c(-1)
));

addGate('G', new Matrix(
    c(1), c(sqrt(2)),
    c(sqrt(2)), c(-1)
).scale(1 / sqrt(3)));

addGate('H', new Matrix(
    c(1), c(1),
    c(1), c(-1)
).scale(1 / sqrt(2)));

addGate('S', new Matrix(
    c(1), c(0),
    c(0), c(0, Math.PI / 2).exp()
));

addGate('T', new Matrix(
    c(1), c(0),
    c(0), c(0, Math.PI / 4).exp()
));

function addGate(name: string, gate: Matrix, tol: number = 1e-5) {
    if(!gate.isUnitary()) {
        console.error('---');
        gate.log();
        throw new Error('addGate: Cannot add non-unitary gate');
    }

    availableGates.push(gate);
    gatesAndNames.push([gate, name]);
    nameToGate[name] = gate;

    const daggered = gate.clone().dagger();

    const daggeredIsSame = daggered.distance(gate) < tol;

    if(!daggeredIsSame) {
        availableGates.push(daggered);
    }

    gatesAndNames.push([daggered, name + '_dagger']);
    nameToGate[name + '_dagger'] = daggered;
}

// Compute all gate set pairs up to n = 5
for(let i = 1; i < 5; i++){
    gateSets.push(...generateGateSets(i));
}

function generateGateSets(n: number = 5): MatrixSet[] {
    if(n === 1) {
        const arr = [];
        for(let i = 0; i < availableGates.length; i++) {
            arr.push(new MatrixSet(availableGates[i]));
        }
        return arr;
    }

    const previousSets = generateGateSets(n - 1);

    const newSets = [];

    for(let i = 0; i < previousSets.length; i++) {
        for(let j = 0; j < availableGates.length; j++) {
            newSets.push(new MatrixSet(previousSets[i], availableGates[j]));
        }
    }

    return newSets;
}

// export const availableGates: Matrix[] = [
//     H,
//     S,
//     S_dagger,
//     T,
//     T_dagger,
//     X,
//     Y,
//     Z,
//     G
// ];

// export const gatesAndNames: [Matrix, string][] = [
//     [H, 'H'],
//     [H, 'H†'],
//     [S, 'S'],
//     [S_dagger, 'S†'],
//     [T, 'T'],
//     [T_dagger, 'T†'],
//     [X, 'X'],
//     [X, 'X†'],
//     [Y, 'Y'],
//     [Y, 'Y†'],
//     [Z, 'Z'],
//     [Z, 'Z†'],
//     [G, 'G'],
//     [G, 'G†'],
// ];

// export const nameToGate: Record<string, Matrix> = {
//     H,
//     S,
//     S_dagger,
//     T,
//     T_dagger,
//     X,
//     Y,
//     Z,
//     G,
//     G_dagger: G
// }