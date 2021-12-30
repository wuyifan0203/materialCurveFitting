const p0 = [5.551383614511336,
  -4.674829278492341,
  0.480299752849433,
  3.372749338864088,
  -1.251536752771733,
  9.134308097402470,
  0.522692455226915,
  -4.240471807705853,
  -28.834770735842053,
  -120.4328257536638,
  2.915192613409641,
  3.518883159679624];

// const h = 6.6260755e-34; // J * s
// const eV_J = 1.6021766208e-19; // 1 eV =1.6021766208e-19 J

// const omega = omega_rad / (2 * Math.PI) * h / eV_J;

const rr1 = p0[0];
const ri1 = p0[1];
const pr1 = p0[2];
const pi1 = p0[3];

const rr2 = p0[4];
const ri2 = p0[5];
const pr2 = p0[6];
const pi2 = p0[7];

const rr3 = p0[8];
const ri3 = p0[9];
const pr3 = p0[10];
const pi3 = p0[11];

const a1 = -(pr1 * pr1);
const b1 = pi1;
const c1 = rr1;
const d1 = ri1;

const a2 = -(pr2 * pr2);
const b2 = pi2;
const c2 = rr2;
const d2 = ri2;

const a3 = -(pr3 * pr3);
const b3 = pi3;
const c3 = rr3;
const d3 = ri3;

const primaryCoe = [
  {
    c: -(pr1 * pr1),
    d: pi1,
    a: rr1,
    b: ri1
  },
  {
    c: -(pr2 * pr2),
    d: pi2,
    a: rr2,
    b: ri2
  },
  {
    c: -(pr3 * pr3),
    d: pi3,
    a: rr3,
    b: ri3
  }
]


module.exports = {
  primaryCoe
}
