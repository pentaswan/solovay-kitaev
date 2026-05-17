import { clampZeroToTwoPi } from "../util";
import { gatesAndNames } from "./availableGates";
import { ComplexNumber } from "./complex";
import { ComplexVector } from "./complexVector";

export class Matrix {
    private data: ComplexNumber[];

    constructor(a: ComplexNumber=new ComplexNumber(1), b: ComplexNumber=new ComplexNumber(0), c: ComplexNumber=new ComplexNumber(0), d: ComplexNumber=new ComplexNumber(1)) {
        this.data = [
            a, b,
            c, d
        ];
    }

    getElement(x: number, y: number): ComplexNumber {
        return this.data[x * 2 + y];
    }

    add(other: Matrix): Matrix {
        for(let i = 0; i < this.data.length; i++) {
            this.data[i].add(other.data[i]);
        }

        return this;
    }

    sub(other: Matrix): Matrix {
        for(let i = 0; i < this.data.length; i++) {
            this.data[i].sub(other.data[i]);
        }

        return this;
    }

    clone(): Matrix {
        return new Matrix(this.data[0].clone(), this.data[1].clone(), this.data[2].clone(), this.data[3].clone());
    }

    trace(): ComplexNumber {
        return this.data[0].clone().add(this.data[2]);
    }

    // computes trace(sqrt(X^t * X)) efficiently
    distance(other: Matrix): number {
        const [a, b, c, d] = this.clone().sub(other).data;

        // |a|^2 + |b|^2 + |c|^2 + |d|^2
        const S = a.magSquared() + b.magSquared() + c.magSquared() + d.magSquared();

        // |ad - bc|
        const D = (a.clone().mul(d).sub(b.clone().mul(c))).mag();

        return Math.sqrt(S + 2 * D);
    }

    dagger(): Matrix {
        const data = this.data;

        this.data = [
            data[0].conjugate(), data[2].conjugate(),
            data[1].conjugate(), data[3].conjugate()
        ]

        return this;
    }

    transpose(): Matrix {
        const data = this.data;

        this.data = [
            data[0], data[2],
            data[1], data[3]
        ]

        return this;
    }

    mul(other: Matrix): Matrix {
        const [a0, a1, a2, a3] = this.data;
        const [b0, b1, b2, b3] = other.data;

        // a0b0 + a1b2, a0b1 + a1b3
        // a2b0 + a3b2, a2b1 + a3b3
        this.data = [
            a0.clone().mul(b0).add(a1.clone().mul(b2)), a0.clone().mul(b1).add(a1.clone().mul(b3)),
            a2.clone().mul(b0).add(a3.clone().mul(b2)), a2.clone().mul(b1).add(a3.clone().mul(b3)),
        ]

        return this;
    }

    scale(scalar: number): Matrix {
        for(let i = 0; i < this.data.length; i++) {
            this.data[i].scale(scalar);
        }

        return this;
    }

    scaleComplex(c: ComplexNumber): Matrix {
        for(let i = 0; i < this.data.length; i++) {
            this.data[i].mul(c);
        }

        return this;
    }

    determinant(): ComplexNumber {
        return (this.data[0].clone().mul(this.data[3]).sub(this.data[1].clone().mul(this.data[2])));
    }

    normalize(): this {
        const det = this.determinant();
        const scalar = 1 / Math.sqrt(det.mag());

        for(let i = 0; i < this.data.length; i++) {
            this.data[i].scale(scalar);
        }

        return this;
    }

    // We want to make the det 1 by multiplying by some complex number c, |c| = 1.
    // C will change the phase twofold because we're multiplying both the top and bottom
    // elements of the matrix by it and the determinant is proportional to the product.
    // So, multiply by e^(-iθ/2) instead of e^(-iθ).
    // Other than that, this function does some tricks if phase is close to pi and that's it.

    transformSU2(tol: number=1e-5): [Matrix, number] {
        const det = this.determinant();
        const detPhase = clampZeroToTwoPi(Math.atan2(det.imaginary(), det.real()));

        const gPhase = detPhase / 2;

        const gPhaseCloseToPi = Math.abs(clampZeroToTwoPi(gPhase) - Math.PI) < tol;
        const rPhase = (-1) ** (gPhaseCloseToPi ? 1 : 0);

        this.scaleComplex(new ComplexNumber(0, -gPhase).exp()).scale(rPhase);

        return [this, rPhase === 1 ? gPhase : 0];
    }

    isUnitary(tol: number=1e-5): boolean {
        const u_dagger_u = this.clone().mul(this.clone().dagger());

        return u_dagger_u.distance(new Matrix()) < tol;
    }

    isHermitian(tol: number=1e-5): boolean {
        // console.log('----');
        // this.log();
        // this.clone().dagger().log();
        // console.log('----');
        return this.distance(this.clone().dagger()) < tol;
    }

    eigenvalues(tol: number=1e-5): [ComplexNumber, ComplexNumber] {
        const [a, b, c, d] = this.data;

        // -(a + d)
        const B = a.clone().add(d).scale(-1);
        // ad - bc
        const C = a.clone().mul(d).sub(b.clone().mul(c));

        // solving λ^2 + Bλ + C = 0
        // sqrt(B^2 - 4AC), but A = 1
        const discriminant = B.clone().mul(B).sub(C.clone().scale(4)).pow(0.5);

        // (-B +- sqrt(B^2 - 4AC)) / 2A
        const λ1 = (B.clone().scale(-1).add(discriminant)).scale(0.5);
        const λ2 = (B.clone().scale(-1).sub(discriminant)).scale(0.5);

        // const arr: [ComplexNumber, ComplexNumber] = [λ1, λ2];

        // arr.sort((a, b) => {
        //     const realDif = a.real() - b.real();
        //     if(Math.abs(realDif) > tol) return Math.sign(realDif);

        //     const imaginaryDif = a.imaginary() - b.imaginary();
        //     if(Math.abs(imaginaryDif) > tol) return Math.sign(imaginaryDif);

        //     return 0;
        // });

        return [λ1, λ2];
    }

    // returns eigenvectors as the columns of a matrix.
    // Assumes the matrix is orthogonal
    eigenvectors(tol: number=1e-5): Matrix {
        if(!this.isUnitary()) {
            console.error('---');
            console.log('Matrix:');
            this.log();
            console.log('Daggered:');
            this.clone().dagger().log();
            console.log('Unitary Product:');
            this.clone().mul(this.clone().dagger()).log();
            throw new Error("Eigenvectors error: Matrix is not unitary!");
        }

        const [λ1, _λ2] = this.eigenvalues();

        const [a, b, _c, _d] = this.data;

        let v1 = new ComplexVector(b.clone(), λ1.clone().sub(a)).normalize();
        // let v2 = v1.clone().orthogonal();

        if(v1.x.mag() < tol) {
            return new Matrix();
        }

        // bigger magnitude should occupy the top left / bottom right 
        if(v1.x.mag() < v1.y.mag()) {
            v1 = v1.orthogonal();
        }

        // rotate by phase into SU(2)
        const phase = Math.atan2(v1.x.imaginary(), v1.x.real());
        const exp = new ComplexNumber(0, -phase).exp();

        v1.x.mul(exp);
        v1.y.mul(exp);

        return new Matrix(
            v1.x, v1.y.clone().conjugate().scale(-1),
            v1.y, v1.x.clone()
        );
    }

    toString(short: boolean = false): string {
        const name = this.name();

        const displayName = name !== '?';
        return `${displayName ? `${name.replace('_dagger', '†')}${short ? '' : ' '}` : ''}${(short && displayName) ? '' : `[
    ${this.data[0].toString()}, ${this.data[1].toString()}
    ${this.data[2].toString()}, ${this.data[3].toString()}
]`}`;
    }

    addNoise(amount: number=0.01, type: NoiseType = NoiseType.Hermitian): Matrix {
        if(type === NoiseType.Hermitian) {
            this.data[0].addNoise(amount).setImaginary(0);
            this.data[1].addNoise(amount);
            this.data[2] = this.data[1].clone().conjugate();
            this.data[3].addNoise(amount).setImaginary(0);
        } else {
            const a = new ComplexNumber(1).addNoise(amount);
            const b = new ComplexNumber(0).addNoise(amount);
            const phi = new ComplexNumber(0).addNoise(amount).setImaginary(0);

            const e_iphi = phi.mul(new ComplexNumber(0, 1)).exp();
            
            // General form for 2x2 unitary matrix:
            // [
            //     a, b,
            //     -e^(iφ)b*, e^(iφ)a*
            // ]
            this.data[0] = a;
            this.data[1] = b;
            this.data[2] = e_iphi.clone().mul(b.clone().conjugate()).scale(-1);
            this.data[3] = e_iphi.clone().mul(a.clone().conjugate());
        }

        this.normalize();

        return this;
    }

    exp(steps: number=15) {
        const X = new Matrix().add(this.clone().scale(1 / (2 ** steps)));
    
        for(let i = 0; i < steps; i++) {
            X.mul(X.clone());
        }

        return X;
    }

    name(tol: number=1e-5): string {
        for(let i = 0; i < gatesAndNames.length; i++) {
            const [gate, name] = gatesAndNames[i];
            if(this.distance(gate) < tol) return name;
        }

        return '?';
    }

    log(short: boolean = false) {
        console.log(this.toString(short));
    }
}

export enum NoiseType {
    Unitary,
    Hermitian
}