const {
    materialParams
} = require('./database')
const {
    primaryCoe
} = require('./p0')

const eps0 = 8.854187817e-12;
const mu0 = 4*Math.PI*1e-7;
var c0 = 1/Math.sqrt(eps0*mu0);
let RwTimes=0;
let XkTimes=0;

    // 极化率函数
const Xk = ({
        w,
        a,
        b,
        c,//g
        d //h
    }) => {
        XkTimes+=1;
        const G2 = c * c;
        const W2 = w * w;;
        const H2 = d * d;
        const ag = a * c;
        const bh = b * d;
        const P = -W2+G2+H2;
        const Y =(-2*(ag+bh));
        const X = 2*a*w;
        const Q = -2*c*w
        const denominator = (-W2+G2+H2)*(-W2+G2+H2)+4*W2*G2
        const real =Y*P+X*Q;
        const virtual = X*P-Y*Q;
        return [real / denominator, virtual / denominator]
}
    // 介电常数
const Rw = (w, Woo,coe) => {
        RwTimes+=1;
        const result = [0, 0]
        for (let i = 0; i < coe.length; i++) {
            const currentCoe = coe[i]
            const XkResult = Xk({
                w,
                a: currentCoe.a,
                b: currentCoe.b,
                c: currentCoe.c,
                d: currentCoe.d
            });
            result[0] += XkResult[0];
            result[1] += XkResult[1]
        }
        result[0] = result[0] + Woo
        return result
}

const materialCurveFitting = ({
    materialId,
    wavelength_min,
    wavelength_max,
    max_of_number = 12,
    weight,
    coe=primaryCoe
}) => {
    const RwList=[]
    let list = materialParams.filter(e => e.wavelength <= wavelength_max && e.wavelength >= wavelength_min);
    const Woo = 1;
    
    let subReal=0;
    let subVirtual=0;

    for (var i=list.length-1;i>=0;i--) {
        let target = list[i];
        // const w = target.wavelength * 2 * Math.PI; //w弧度单位
        const currentWavelength=target.wavelength*1000000 //转换微米
        const w = 2*Math.PI*c0/ (currentWavelength * 1e-6)
        const  h = 6.6260755e-34; // J * s
        const eV_J = 1.6021766208e-19; // 1 eV =1.6021766208e-19 J

        const omega = w / (2*Math.PI) * h /eV_J;//电子伏单位
        // omega_rad = 2 * pi * c0 ./ (lambda * 1e-6);

        const RwResult = Rw(omega,Woo,coe);
        //折射率
        const ReIndex = 
        RwList.push(RwResult);
        // console.log(RwResult);
        // console.log([target.wavelength*1000000,target.epsIm]);
        //   console.log([target.wavelength*1000000,RwResult[1]]);
       
        const DetaReal = RwResult[0]-target.epsRe;
        const DetaVirtual = RwResult[1]+target.epsIm;
        // console.log(DetaVirtual);
        subReal+=DetaReal*DetaReal;
        subVirtual+= DetaVirtual*DetaVirtual;
    }
    // console.log('---------------------start------------------------------');
    // console.log("subreal",subReal);
    // console.log("subVirtual",subVirtual);

    const RMSE = Math.sqrt((subVirtual+subReal)/list.length);
    const weightedRMSE =Math.sqrt((weight*subVirtual+subReal)/list.length);
    // console.log('list length:',list.length);
    // console.log('Xk Times:',XkTimes);
    // console.log('Rw Times:',RwTimes);
    // console.log('RMSE :',RMSE);
    // console.log('-----------------------end------------------------------');
    return {RMSE,weightedRMSE,RwList}
};

let rmse=materialCurveFitting({
    wavelength_min: 0.2e-6,
    wavelength_max: 2.0e-6,
    weight: 1,
}).RMSE
// console.log(rmse);


module.exports = {materialCurveFitting};