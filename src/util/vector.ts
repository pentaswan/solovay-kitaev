export class Vector {
    public data: number[];
    constructor(x: number = 0, y: number = 0, z: number = 0) {
        this.data = [x, y, z];
    }

    public get x() {
        return this.data[0];
    }

    public get y() {
        return this.data[1];
    }

    public get z() {
        return this.data[2];
    }

    toString() {
        return `(${this.data.join(', ')})`;
    }

    log() {
        console.log(this.toString());
    }
}