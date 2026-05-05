import { toMinimumPrecision } from "../util";

export class ComplexNumber {
    private x: number;
    private y: number;

    constructor(x: number = 0, y: number = 0) {
        this.x = x;
        this.y = y;
    }

    set(x: number = 0, y: number = 0) {
        this.x = x;
        this.y = y;
    }

    add(other: ComplexNumber): ComplexNumber {
        this.x += other.x;
        this.y += other.y;

        return this;
    }

    sub(other: ComplexNumber): ComplexNumber {
        this.x -= other.x;
        this.y -= other.y;

        return this;
    }

    real(): number {
        return this.x;
    }

    imaginary(): number {
        return this.y;
    }

    setReal(num: number): ComplexNumber {
        this.x = num;

        return this;
    }

    setImaginary(num: number): ComplexNumber {
        this.y = num;

        return this;
    }

    mul(other: ComplexNumber): ComplexNumber {
        // (x + iy) * (other.x + i*other.y)
        const newX = this.x * other.x - this.y * other.y;
        const newY = this.y * other.x + this.x * other.y;

        this.x = newX;
        this.y = newY;

        return this;
    }

    mag(): number {
        return Math.sqrt(this.x ** 2 + this.y ** 2);
    }

    magSquared(): number {
        return this.x ** 2 + this.y ** 2;
    }

    scale(scalar: number): ComplexNumber {
        this.x *= scalar;
        this.y *= scalar;

        return this;
    }

    conjugate(): this {
        this.y *= -1;

        return this;
    }

    interpolate(other: ComplexNumber, t: number): ComplexNumber {
        this.x = this.x * (1 - t) + other.x * t;
        this.y = this.y * (1 - t) + other.y * t;

        return this;
    }

    clone(): ComplexNumber {
        return new ComplexNumber(this.x, this.y);
    }

    exp(): ComplexNumber {
        const exp = Math.exp(this.x);

        const newX = exp * Math.cos(this.y);
        const newY = exp * Math.sin(this.y);

        this.x = newX;
        this.y = newY;

        return this;
    }

    pow(power: number) {
        const r = this.mag();
        const theta = Math.atan2(this.y, this.x);

        const newR = r ** power;
        const newTheta = theta * power;

        this.x = newR * Math.cos(newTheta);
        this.y = newR * Math.sin(newTheta);

        return this;
    }

    addNoise(amount: number=0.01): ComplexNumber {
        this.x += (Math.random() * 2 - 1) * amount;
        this.y += (Math.random() * 2 - 1) * amount;

        return this;
    }

    toString(tol: number=1e-5): string {
        let x = toMinimumPrecision(this.x, 3);
        let y = toMinimumPrecision(Math.abs(this.y), 3);

        let xIsZero = Math.abs(this.x) < tol;
        let yIsZero = Math.abs(this.y) < tol;
        let yIsOne = Math.abs(Math.abs(this.y) - 1) < tol;

        if(xIsZero && yIsZero) {
            return '0';
        } 
        
        // x = 0, y != 0. Form a * i
        if(xIsZero) {
            return `${this.y >= 0 ? '' : '-'}${yIsOne ? '' : y}i`;
        }

        // x != 0, y = 0. Form a.
        if(yIsZero) {
            return x;
        }

        // x != 0, y != 0
        return `${x}${this.y >= 0 ? '+' : '-'}${y}i`;
    }

    log() {
        console.log(this.toString());
    }
}