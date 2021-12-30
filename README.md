<a name="LvXmQ"></a>
# 需要了解的公式
<a name="E7r1J"></a>
## 公式一
真空介电常量是物理量在度量时引进的常数(主要是库仑定律中对电荷量的度量)，根据麦克斯韦方程组，可推知真空介电常数与其它物理常数的关系：<br />![](https://cdn.nlark.com/yuque/0/2021/jpeg/25638087/1640829257428-35527a13-a754-4911-8e50-6b7cf8037a25.jpeg#clientId=uf9577f2a-5334-4&crop=0&crop=0&crop=1&crop=1&from=paste&id=u543553da&margin=%5Bobject%20Object%5D&originHeight=46&originWidth=81&originalType=url&ratio=1&rotation=0&showTitle=false&status=done&style=none&taskId=u74694dac-5428-4da7-8b29-c0759647aae&title=)<br />其中，<br />![](https://cdn.nlark.com/yuque/0/2021/jpeg/25638087/1640829359762-a9a53ddd-d129-4f1e-9cf4-aa5360f9362e.jpeg#clientId=uf9577f2a-5334-4&crop=0&crop=0&crop=1&crop=1&from=paste&id=u7dc61bac&margin=%5Bobject%20Object%5D&originHeight=21&originWidth=157&originalType=url&ratio=1&rotation=0&showTitle=false&status=done&style=none&taskId=u26d3cfd2-35d4-4cfa-869d-104c1b835b5&title=)<br />_ε0为真空介电常数，值为8.854187817*10^-12 F/m_<br />c为光速<br />

<a name="a5fNk"></a>
## 公式二
光速与波长，频率的关系<br />![image.png](https://cdn.nlark.com/yuque/0/2021/png/25638087/1640852838704-75cac76f-c686-456f-80d3-699881e0438b.png#clientId=u3c323c1c-af86-4&crop=0&crop=0&crop=1&crop=1&from=paste&height=48&id=u2c02f86c&margin=%5Bobject%20Object%5D&name=image.png&originHeight=48&originWidth=81&originalType=binary&ratio=1&rotation=0&showTitle=false&size=1011&status=done&style=none&taskId=ue2c6537f-4a24-43eb-86e6-db71b2d7dd2&title=&width=81)<br />其中，c为光速，f为频率，λ为波长。<br />​<br />
<a name="WJeK6"></a>
## 公式三
E = hv<br />其中，E为能量，h为普朗克常数，v为频率。<br />h = 6.6260755e-34 J·s<br />​<br />
<a name="bOkda"></a>
## 公式四
RMSE[均方根误差公式](https://baike.baidu.com/item/%E5%9D%87%E6%96%B9%E6%A0%B9%E8%AF%AF%E5%B7%AE/3498959?fromtitle=RMSE&fromid=6536667&fr=aladdin)<br />RMSE = $$\sqrt{\frac{1}{m}\sum_{n=1}^m |y_d - y_i|^2 }$$<br />​<br />
<a name="BGAxE"></a>
## 公式五
介电常数与极化率公式<br />![企业微信截图_16400837131913.png](https://cdn.nlark.com/yuque/0/2021/png/25638087/1640831678381-fa15bdd2-2299-49e5-b070-f8b48c847502.png#clientId=uf9577f2a-5334-4&crop=0&crop=0&crop=1&crop=1&from=paste&height=137&id=u0c27af26&margin=%5Bobject%20Object%5D&name=%E4%BC%81%E4%B8%9A%E5%BE%AE%E4%BF%A1%E6%88%AA%E5%9B%BE_16400837131913.png&originHeight=308&originWidth=830&originalType=binary&ratio=1&rotation=0&showTitle=false&size=42058&status=done&style=none&taskId=u5329c961-de2a-4b4e-9a88-77ddb2a7f03&title=&width=368)<br />
<br />接下来，就可以正式开始材料拟合了<br />​<br />
<a name="zHMDY"></a>
# 代码书写
注：本次拟合材料以Si为例
<a name="qdJNv"></a>
## 步骤
将整体分为三步🚶‍♀️计算

1. 计算极化率（公式五的Xk）；
1. 计算拟合点（公式一的_εr_）；
1. 计算均方根误差；



<a name="lDPRu"></a>
## 代码
<a name="PrNi0"></a>
### 1.计算极化率
![image.png](https://cdn.nlark.com/yuque/0/2021/png/25638087/1640853334198-57735fb2-2d6b-4870-aa0d-9c30aeabc85e.png#clientId=uec2524fb-2b0b-4&crop=0&crop=0&crop=1&crop=1&from=paste&height=48&id=u1e775c2c&margin=%5Bobject%20Object%5D&name=image.png&originHeight=54&originWidth=68&originalType=binary&ratio=1&rotation=0&showTitle=false&size=2256&status=done&style=none&taskId=ub6d9e403-2b99-4409-b427-415f2e9409b&title=&width=60)![image.png](https://cdn.nlark.com/yuque/0/2021/png/25638087/1640846354672-12714f39-62d5-4a1a-981f-9e53558676e8.png#clientId=ua0a95876-f94c-4&crop=0&crop=0&crop=1&crop=1&from=paste&height=45&id=WgqfR&margin=%5Bobject%20Object%5D&name=image.png&originHeight=45&originWidth=137&originalType=binary&ratio=1&rotation=0&showTitle=false&size=5270&status=done&style=none&taskId=u3f2cad2b-4b5b-4bc6-bda6-4e895e0cf8d&title=&width=137)<br />rk、rk*、pk、pk*为复数，rk与rk*，pk与pk*互为共轭复数<br />由于JavaScript没有复数类型，进行拆解，方便计算<br />设

- rk  = a + bi
- rk* = a - bi
- pk = c + di
- pk* = c - di

接下来进行通分，比较复杂恶心，且容易出错，就不在这里化简了<br />​

最后化简好，将实部虚部化成形如格式 为`X + Y i `就可以了<br />入参为(a,b,c,d,w) 返回值为数组分别存储 `X`与`Y`，<br />第一步就完成了，代码如下；
```javascript
const Xk = ({w,a,b,c,d }) => {
        const [C2,W2,D2] = [c * c,w * w,d * d];
        const ac = a * c;
        const bd = b * d;
        const P = -W2+C2+H2;
        const Y =(-2*(ac+bd));
        const X = 2*a*w;
        const Q = -2*c*w
        const denominator = (-W2+C2+D2)*(-W2+C2+D2)+4*W2*C2;
        const real =Y*P+X*Q;
        const imag = X*P-Y*Q;
        return [real / denominator, imag / denominator]
}
```


<a name="Mtmr2"></a>
### 2.计算拟合点
![image.png](https://cdn.nlark.com/yuque/0/2021/png/25638087/1640848417739-fc7e51d0-59e0-4b4b-aaf3-3dd88f4e497d.png#clientId=ua0a95876-f94c-4&crop=0&crop=0&crop=1&crop=1&from=paste&height=54&id=ud4ab786f&margin=%5Bobject%20Object%5D&name=image.png&originHeight=54&originWidth=190&originalType=binary&ratio=1&rotation=0&showTitle=false&size=6642&status=done&style=none&taskId=u9069b2f6-a459-4a0c-8947-2dd6e96b521&title=&width=190)<br />![image.png](https://cdn.nlark.com/yuque/0/2021/png/25638087/1640848319800-832ce820-e56c-4f5c-98ef-0d4ae4b05003.png#clientId=ua0a95876-f94c-4&crop=0&crop=0&crop=1&crop=1&from=paste&height=14&id=ua6be757e&margin=%5Bobject%20Object%5D&name=image.png&originHeight=14&originWidth=22&originalType=binary&ratio=1&rotation=0&showTitle=false&size=698&status=done&style=none&taskId=u2e47229a-a783-4066-becf-f7c59d2eaa9&title=&width=22)=1<br />默认三对极点 L=3<br />将三组系数存储在coe中
```javascript
const coe = [
  {
    a:xxxx,
    b:xxxx,
    c:xxxx,
    d:xxxx
  },
  ...
]
```
入参为(coe,Woo,w) 返回值依然为复数数组分别存储，<br />第二步就完成了，代码如下：
```javascript
const Rw = (w, Woo,coe) => {
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
```


<a name="v1YrW"></a>
### 3.计算均方根误差
RMSE = $$\sqrt{\frac{1}{m}\sum_{n=1}^m |y_d - y_i|^2 }$$<br />m为点的个数<br />通过第二步，已经求出了拟合点，这时候需要计算一个点的误差<br />`单点误差 = 实际测量点 - 拟合点`<br />但由于实际测量点与拟合点都是复数，所以采用<br />`单点实部误差 = 实际测量点实部 - 拟合点实部`<br />`单点虚部误差 = 实际测量点虚部 - 拟合点虚部`<br />最后在平方求和，除以点数，开根号，即可完成计算<br />​

计算加权误差<br />`加权误差 = 误差虚部的平方的和 x 权重`<br />​

入参为({materialId,wavelength_min,wavelength_max,weight,coe,max_of_number}) <br />返回值为对象{RMSE,weightedRMSE,RwList}，<br />第三步就完成了，代码如下：
```javascript
const materialCurveFitting = ({
    materialId,
    wavelength_min,
    wavelength_min,
    max_of_number = 6, // 极点个数，三对极点
    weight,
    coe=primaryCoe
}) => {
    //从数据库提取所需材料的符合wavelength_min，wavelength_max之间所有的点
    // 注意 ⚠️ 应该为闭区间
    let list = materialParams.filter(e => e.wavelength <= wavelength_max && e.wavelength >= wavelength_min);
    //
    const Woo = 1;
    
    let subReal=0; //实部总和初始化
    let subVirtual=0; //虚部总和初始化
    const RwList=[]；//所有拟合点的实部虚部的集合

    for (var i=list.length-1;i>=0;i--) {
        let target = list[i];
        const currentWavelength=target.wavelength*1000000 //转换微米
        //公式一
        const w = 2*Math.PI*c0/ (currentWavelength * 1e-6)
        const h = 6.6260755e-34; // J * s
        //公式三
        const eV_J = 1.6021766208e-19; // 1 eV =1.6021766208e-19 
        //公式二
        const omega = w / (2*Math.PI) * h /eV_J;//电子伏单位
        
        //计算拟合点
        const RwResult = Rw(omega,Woo,coe);
        RwList.push(RwResult);
       
        const DetaReal = RwResult[0]-target.epsRe;
        const DetaVirtual = RwResult[1]+target.epsIm;
        //求和
        subReal+=DetaReal*DetaReal;
        subVirtual+= DetaVirtual*DetaVirtual;
    }
    // console.log('---------------------start------------------------------');
    // console.log("subreal",subReal);
    // console.log("subVirtual",subVirtual);
  
    // RMSE 的计算
    const RMSE = Math.sqrt((subVirtual+subReal)/list.length);
    // 加权后的RMSE
    const weightedRMSE =Math.sqrt((weight*subVirtual+subReal)/list.length);
    
    // console.log('list length:',list.length);
    // console.log('RMSE :',RMSE);
    // console.log('-----------------------end------------------------------');
    return {RMSE,weightedRMSE,RwList}
};
```
其中c0通过公式一计算
```javascript
const eps0 = 8.854187817e-12;
const mu0 = 4*Math.PI*1e-7;
var c0 = 1/Math.sqrt(eps0*mu0);
```


<a name="ZQtdy"></a>
# 函数调用
```javascript
const { RMSE, weightedRMSE, RwList } = materialCurveFitting({
    wavelength_min: 0.2e-6,
    wavelength_max: 2.0e-6,
    weight: 1,
})
console.log('拟合误差：',RMSE);
console.log('加权拟合误差：',weightedRMSE);
console.log('拟合点集合：',RwList);
```
​<br />
<a name="uonmH"></a>
# 致谢
非常感谢技术部 _刘琦 _老师提供技术指导以及文档的订正

