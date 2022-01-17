//Nelder_Mead.m
const {materialCurveFitting} = require("./demo.js")
const mathjs = require('mathjs')

const Tol = 0.1;// 目标误差
const MaxIter = 10000
;

const START = Date.now();

let counter = 0;
let err = 10000; //误差

// const p0 = [
//   5.551383614511336,
//   -4.674829278492341,
//   0.480299752849433,
//   3.372749338864088,
//   -1.251536752771733,
//   9.134308097402470,
//   0.522692455226915,
//   -4.240471807705853,
//   -28.834770735842053,
//   -120.4328257536638,
//   2.915192613409641,
//   3.518883159679624
// ];

const p0 =[
  5.700686837863513,
4.578720280401709,
0.482102135864590,
-3.369776157658570,
-1.317293027611639,
-8.835480500791704,
0.516962791801114,
4.243456713267605,
-17.426872989705963,
-50.139855115730576,
2.430251701626421,
4.596657303636679
]

// const makeP =()=>{
//   const p=[];
//   for(let i=0;i<12;i++){
//     p.push(Math.floor(Math.random() * 20+ 1));
//   }
//   return p
// }
// const p0=makeP();
// console.log(p0);

// const p0=[
//   6, 14, 12, 10, 14,
//   2,  9,  5,  1, 14,
//   8, 15
// ]

// const p0 = [-1.697748954791354,
//       2.219730290176961e+03,
//       0.238345515816527,
//       -0.017013360108035,
//       1.347382972425653,
//       -0.762209952188253,
//       -0.649660513177193,
//       2.584944231946596,
//       -4.417405293428975e+04,
//       -1.910111189502858e+06,
//       16.429865129839193,
//       6.338636303438941];

const lambada = [
  {
    rr:1.1e1,
    ri:1.2e1,
    pr:1.3e1,
    pi:1.4e1
  },
  {
    rr:1.5e1,
    ri:1.6e1,
    pr:1.7e1,
    pi:1.8e1
  },
  {
    rr:1.9e1,
    ri:2.0e1,
    pr:2.1e1,
    pi:2.2e1
  }
];

const optimizationLengthScale =[
  1.1e1,
  1.2e1,
  1.3e1,
  1.4e1,
  1.5e1,
  1.6e1,
  1.7e1,
  1.8e1,
  1.9e1,
  2.0e1,
  2.1e1,
  2.2e1
]
// 二维数组创建初始化
const zeros = (n, m) => {
  const tArray = new Array(); //先声明一维
  for (let i = 0; i < n; i++) { //一维长度为i,i为变量，可以根据实际情况改变
    tArray[i] = new Array(); //声明二维，每一个一维数组里面的一个元素都是一个数组；
    for (let j = 0; j < m; j++) { //一维数组里面每个元素数组可以包含的数量p，p也是一个变量；
      tArray[i][j] = 0; //这里将变量初始化，我这边统一初始化为空，后面在用所需的值覆盖里面的值
    }
  }
  return tArray
}

const subMatrix = (Matrix) =>{
  let result= mathjs.zeros(Matrix[0].length)
  for(let i=0;i<Matrix.length;i++){
    result = mathjs.add(result,Matrix[i])
  }
  return result;
}

//生成对角矩阵
const diag =  (Array) => {
  const p = zeros(Array.length,Array.length);
  for(let i=0;i<Array.length;i++){
    p[i][i]=Array[i]
  }
  return p
}
// 生成一组
const pkrk = ({rr,ri,pr,pi})=>{
  return {
    a:rr,
    b:ri,
    c:-(pr*pr),
    d:pi
  }
}
// 生成一行三组
const pkrkColumn = (Plist)=>{
  let j=0
  const pkrkcolumn=[]
  for(let i=0;i<Plist.length/4;i++){
    const p = pkrk({
      rr:Plist[j*4],
      ri:Plist[j*4+1],
      pr:Plist[j*4+2],
      pi:Plist[j*4+3]
    })
    j=j+1;
    pkrkcolumn.push(p)
  }
  return pkrkcolumn;
}
// 更新Plist
const updatePMatrix =(indexArray,PMatrix) =>{
  const newPmatrix = [];
  for(let i=0;i<indexArray.length;i++){
    const index=indexArray[i].index;
    newPmatrix[i]=PMatrix[index];
  }
  return newPmatrix
}
//生成初始化P矩阵
const initPMatrix=  (p0,optimizationLengthScale)=>{
  const p=  diag(optimizationLengthScale);
  for(let i=0;i<p.length;i++){
    for(let j=0;j<p0.length;j++){
      p[i][j]+=p0[j]
    }
  }
  p.unshift(p0);
  return p
};

// function quick(arr){
//   if(arr.length<=1){
//       return arr;
//   }
//   let left = [],
//       right = [],
//       k = arr.splice(0, 1)[0].value;
//       for(var i = 0;i<arr.length;i++){
//           if(arr[i].value>k){
//               right.push(arr[i].value);
//           }else{
//               left.push(arr[i].value);
//           }
//       }
//       return quick(left).concat([k],quick(right));
// }

const sortArray=(arr)=>{
  const A=[];
  for(let i=0;i<arr.length;i++){
    A.push({value:arr[i],index:i})
  }
  if(A.length<1){
    return
  }
    let len = A.length;
    
    for (let i = 0; i < len-1; i++) {
      for (let j = 0; j < len - 1 - i; j++) {
          if (A[j].value > A[j + 1].value) {
              let temp = A[j];
              A[j] = A[j+1];
              A[j+1] = temp;
          }
      }
    }
    return A;
};

const arrayMultiplier=(num,array)=>{
  for(let i=0;i<array.length;i++){
    array[i] = array[i]*num;
  }
  return array;
}

let PMatrix=initPMatrix(p0,optimizationLengthScale);
// const measure = [1e-50,50.0e-6]
const measure = [2e-7,2e-6];
const err_cf = zeros(PMatrix.length,1);
let p00 = PMatrix[0];
const RMSEList = [];
let p_c

while (counter < MaxIter && err > Tol){
  for (let i = 0; i < PMatrix.length; i++) {
    err_cf[i] = materialCurveFitting({
      coe:pkrkColumn(PMatrix[i]),
      wavelength_min:measure[0],
      wavelength_max:measure[1],
    }).RMSE;
  }


  
  let f_p = sortArray(err_cf);
  // let oldP = PMatrix;

  let P_tem =  updatePMatrix(f_p,PMatrix);

  PMatrix = P_tem;



  // console.log(PMatrix);

  let p_s = PMatrix[0];
  let p_l = PMatrix[PMatrix.length-1];
  // console.log('pl:',p_l);
  // let p_nl =PMatrix[PMatrix.length-2];


  let f_s = f_p[0].value;
  let f_l = f_p[f_p.length-1].value;
  let f_nl = f_p[f_p.length-2].value;

  // console.log("fs",f_s);
  // console.log("fl",f_l);
  // console.log("fnl",f_nl);

  // stage 2
  let rho = 1;
  let x = 2;
  let r = 0.5;
  let sigma = 0.5;
  // console.log("PMatrix:",PMatrix);
  // console.log("length",PMatrix.length-1);

  let p_g=arrayMultiplier(1/(PMatrix.length-1),mathjs.subtract(subMatrix(PMatrix)._data,p_l));
  // console.log("subP:",mathjs.subtract(subMatrix(PMatrix)._data,p_l));
  // console.log("p_g:",p_g);
  let p_r=mathjs.add(p_g,arrayMultiplier(rho,mathjs.subtract(p_g,p_l)));
  // console.log("pr",p_r);
  // console.log('复数:',mathjs.sqrt(mathjs.complex(2, 3) ));
  // console.log('nOneMatrix:',OneNMatrix(1/12,12));

  // stage 2.1

  // console.log(p_r,p_l);

  let f_p_r =  materialCurveFitting({
    coe:pkrkColumn(p_r),
    wavelength_min:measure[0],
    wavelength_max:measure[1],
  }).RMSE;

  // console.log("f_p_r:",f_p_r);


  if(f_p_r>=f_s && f_p_r<f_nl){
    // console.log(1);
    PMatrix[PMatrix.length-1]=p_r;

  }else if(f_p_r < f_s){
    // console.log(2);
    let p_e = mathjs.add(p_g,arrayMultiplier(x,mathjs.subtract(p_r, p_g)));
    // console.log(p_e);
    let f_p_e =  materialCurveFitting({
      coe:pkrkColumn(p_e),
      wavelength_min:measure[0],
      wavelength_max:measure[1],
    }).RMSE;
    // console.log('fpe:',f_p_e);
    if(f_p_e < f_p_r){
      // console.log(2.1);
      PMatrix[PMatrix.length-1]=p_e;
    }
    // else if(f_p_e >= f_p_r){
      else{
      // console.log(2.2);
      PMatrix[PMatrix.length-1]=p_r;
    }
  }else if(f_p_r >= f_nl){
    // console.log(3);
    if (f_p_r < f_l){
      // console.log(3.1);
      p_c = mathjs.add(p_g,arrayMultiplier(r,mathjs.subtract(p_r,p_g)));
    }
    // else if(f_p_r >= f_l){
      else{
      // console.log(3.2);
      p_c = mathjs.add(p_g,arrayMultiplier(r,mathjs.subtract(p_l,p_g))); 
    }
    let f_c= materialCurveFitting({
      coe:pkrkColumn(p_c),
      wavelength_min:measure[0],
      wavelength_max:measure[1],
    }).RMSE;
    if (f_c <= f_l){
      // console.log(3.3);
      PMatrix[PMatrix.length-1]=p_c;
    }
    else{
    // else if(f_c > f_l){
      // console.log(3.4);
      for(let i=1;i<PMatrix.length;i++){
        PMatrix[i]=mathjs.add(p_s,arrayMultiplier(sigma,mathjs.subtract(PMatrix[i],p_s)));
      }
    }
  }

  let Rw2List=[];
   let RwList=  materialCurveFitting({
      coe:pkrkColumn(PMatrix[0]),
      wavelength_min:measure[0],
      wavelength_max:measure[1],
   }).RwList;
   for(let i=0;i<RwList.length;i++){
     Rw2List.push(-1*RwList[i][1]);
   }
  Rw2List=Rw2List.sort();
  // console.log('low',Rw2List[0]);
  // console.log('lastPMatrix:',PMatrix[0]);

  if(Rw2List[0]>=0){
    // console.log(4);
    p00=PMatrix[0]
    RMSEList[counter]=  materialCurveFitting({
      coe:pkrkColumn(p00),
      wavelength_min:measure[0],
      wavelength_max:measure[1],
    }).RMSE;
  }else{
    // console.log(5);
    if (counter > 1){
      RMSEList[counter] = RMSEList[counter-1];
    }else{
      RMSEList[counter] = f_s;
    }
  }
  counter++;
}
const END = Date.now();

const SPEND = END-START
console.log(SPEND/1000 +" s");

console.log(RMSEList[RMSEList.length-1]);
// console.log(p00);
// console.log(err);
 