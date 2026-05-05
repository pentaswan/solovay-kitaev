import { ComplexNumber } from "./complex";
import { Matrix } from "./matrix";

export class MatrixSet {
    public matricies: Matrix[];
    public computed: Matrix;

    constructor(...matricies: (Matrix | MatrixSet)[]) {
        this.matricies = [];
        this.computed = new Matrix();

        for(let i = 0; i < matricies.length; i++) {
            const m = matricies[i];
            if(m instanceof MatrixSet) {
                this.matricies.push(...m.matricies.map(m => m.clone()));
                this.computed.mul(m.computed);
            } else if(m instanceof Matrix){
                this.matricies.push(m.clone());
                this.computed.mul(m);
            }
        }
    }

    mul(other: MatrixSet): MatrixSet {
        this.checkComputed();
        
        this.matricies.push(...other.matricies.map(m => m.clone()));

        for(let i = 0; i < other.matricies.length; i++) {
            this.computed.mul(other.matricies[i]);
        }

        this.checkComputed();

        return this;
    }

    clone() {
        const newSet = new MatrixSet();

        newSet.computed = this.computed.clone();
        
        for(let i = 0; i < this.matricies.length; i++) {
            newSet.matricies.push(this.matricies[i].clone());
        }

        return newSet;
    }

    dagger(): MatrixSet {
        this.checkComputed();

        this.computed.dagger();

        // (AB)† = B†A†
        const newMatricies = [];
        for(let i = this.matricies.length - 1; i >= 0; i--) {
            newMatricies.push(this.matricies[i].dagger());
        }

        this.matricies = newMatricies;

        this.checkComputed();

        return this;
    }

    reverse(): MatrixSet {
        const newMatricies = [];
        for(let i = this.matricies.length - 1; i >= 0; i--) {
            newMatricies.push(this.matricies[i]);
        }

        // redo computed, there's no easy way to know it via dagger.
        const ret = new MatrixSet(...newMatricies);

        return ret;
    }

    transpose(): MatrixSet {
        const newMatricies = [];
        for(let i = this.matricies.length - 1; i >= 0; i--) {
            newMatricies.push(this.matricies[i].transpose());
        }

        // redo computed, there's no easy way to know it via dagger.
        const ret = new MatrixSet(...newMatricies);

        return ret;
    }

    mulPhase(phase: number): MatrixSet {
        const exp = new ComplexNumber(0, phase).clone().exp();
        const phaseMat = new MatrixSet(new Matrix(
            exp.clone(), new ComplexNumber(),
            new ComplexNumber(), exp.clone()
        ));

        this.mul(phaseMat);

        return this;
    }

    eigenvectors(): MatrixSet {
        return new MatrixSet(this.computed.eigenvectors());
    }

    log(short: boolean = true) {
        console.log('--- MatrixSet Log ---');

        let str = '';
        for(let i = 0; i < this.matricies.length; i++) {
            const newStr = this.matricies[i].toString(short).replaceAll('\n', '').replace('    ','').replace('    ', ', ');

            let prefix = '';
            if(i !== 0) prefix = ', ';

            str += prefix + newStr;
        }
        console.log(str);
        console.log('Computed:');
        this.computed.log(short);
        console.log('--- End Log ---');
    }

    checkComputed(tol: number = 1e-5): boolean {
        const m = new Matrix();
        for(let i = 0; i < this.matricies.length; i++) {
            m.mul(this.matricies[i]);
        }

        const distance = m.distance(this.computed);

        const valid = distance < tol;

        if(!valid) {
            console.error('---');
            this.log();
            console.log('computed matrix should be:');
            m.log();
            throw new Error("checkComputed error: computed matrix does not match stored matricies.");
        }

        // console.log('👍');

        return valid;
    }
}