export function roundTo(num: number, to: number = 100) {
    return Math.round(num / to) * to;
}

export function clamp(val: number, min: number, max: number) {
    return Math.max(min, Math.min(max, val));
}

export function toMinimumPrecision(num: number, digits: number = 3) {
    let str = num.toPrecision(digits);

    while (str.length > 1) {
        const lastChar = str[str.length - 1];

        if (lastChar !== '0' && lastChar !== '.') break;

        str = str.slice(0, str.length - 1);
    }

    return str;
}

export function clampZeroToTwoPi(angle: number) {
    return (angle + Math.PI * 2) % (Math.PI * 2);
}