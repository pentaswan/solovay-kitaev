import { ComplexNumber } from "./complex";

export class ComplexVector {
    public data: ComplexNumber[];
    constructor(x: ComplexNumber = new ComplexNumber(), y: ComplexNumber = new ComplexNumber(), ...otherCoords: ComplexNumber[]) {
        this.data = [x, y];
        for(let i = 0; i < otherCoords.length; i++) {
            this.data.push(otherCoords[i]);
        }
    }

    public get x(): ComplexNumber {
        return this.data[0];
    }

    public get y(): ComplexNumber {
        return this.data[1];
    }

    public get z(): ComplexNumber {
        return this.data[2];
    }
    
    clone() {
        return new ComplexVector(...this.data.map(a => a.clone()));
    }

    mag(): number {
        let mag = 0;

        for(let i = 0; i < this.data.length; i++) {
            mag += this.data[i].magSquared();
        }

        return Math.sqrt(mag);
    }

    orthogonal(): ComplexVector {
        if(this.data.length !== 2) {
            throw new Error("Orthogonal error: Can only make complex vectors of length 2 orthogonal.");
        }
        
        this.data = [
            this.data[1].scale(-1).conjugate(),
            this.data[0].conjugate()
        ]

        return this;
    }

    normalize(): ComplexVector {
        const mag = this.mag();

        if(mag === 0) return this;

        for(let i = 0; i < this.data.length; i++) {
            this.data[i].scale(1 / mag);
        }

        return this;
    }

    simplifyPhase(tol: number=1e-5): ComplexVector {
        // we can get things like (0, e^(iθ)) out of our eigenvector calculation.
        // Eigenvectors are indifferent to scalar multiplication, so we'll multiply by e^(-iθ) to get simplified eigenvalues.
        let phase = 0;

        if(this.x.mag() > tol) phase = Math.atan2(this.x.imaginary(), this.x.real());
        else if(this.y.mag() > tol) phase = Math.atan2(this.y.imaginary(), this.y.real());

        const phaseShift = new ComplexNumber(0, -phase).exp();

        this.x.mul(phaseShift);
        this.y.mul(phaseShift);

        return this;
    }
}