(function(root){
const rad=a=>a*Math.PI/180, dot=(a,b)=>a.x*b.x+a.y*b.y, sub=(a,b)=>({x:a.x-b.x,y:a.y-b.y}), cross=(a,b)=>a.x*b.y-a.y*b.x;
// World lengths are millimetres (1000 world units = 1 metre).
const isLens=o=>o.type==='lens'||o.type==='diverging';
function lensPower(n,r,t,negative=false){if(!Number.isFinite(r))return 0;return 1000*(n-1)*((negative?-2:2)/r-(n-1)*t/(n*r*r));}
function lensGeometry(o){
 const a=o.h/2,edge=8,negative=o.power<0,minR=Math.max(a*1.05,8),n=o.n;
 if(!Number.isFinite(o.power)||!Number.isFinite(n)||n<1||n>3||o.h<8||o.h>1200)throw Error('Некорректные параметры линзы.');
 const sag=r=>a*a/(r+Math.sqrt(r*r-a*a));
 const powerAt=r=>lensPower(n,r,negative?edge:edge+2*sag(r),negative);
 const max=Math.abs(powerAt(minR));
 if(o.power!==0&&Math.abs(o.power)<0.001)throw Error('Минимальная ненулевая сила — 0,001 дптр. Для пластинки введите 0.');
 if(Math.abs(o.power)>max+1e-8)throw Error('Для этой высоты и n допустима сила до '+max.toFixed(2)+' дптр по модулю. Уменьшите силу или высоту.');
 if(o.power===0)return {radius:Infinity,w:edge,t:edge,max};
 let lo=minR,hi=Math.max(minR*2,4000*(n-1)/Math.abs(o.power));
 for(let i=0;i<80;i++){let mid=(lo+hi)/2;if(Math.abs(powerAt(mid))>Math.abs(o.power))lo=mid;else hi=mid;}
 let radius=(lo+hi)/2,w=edge+2*sag(radius);return {radius,w,t:negative?edge:w,max};
}
function normalizeLens(o){if(!isLens(o))return {...o};let result={...o};
 if(result.power===undefined){let a=o.h/2,sag=Math.min(Math.max((o.w-8)/2,0),a*.7),r=sag>0?(a*a+sag*sag)/(2*sag):Infinity;result.power=lensPower(o.n,r,o.type==='diverging'?8:8+2*sag,o.type==='diverging');}
 result.power=(result.type==='diverging'?-1:1)*Math.abs(result.power);const g=lensGeometry(result);result.w=g.w;return result;
}
function polygon(o){let p=[],w=o.w/2,h=o.h/2;
if(o.type==='circle')for(let i=0;i<160;i++){let a=i*2*Math.PI/160;p.push({x:w*Math.cos(a),y:h*Math.sin(a)});}
else if(o.type==='triangle')p=[{x:0,y:-h},{x:w,y:h},{x:-w,y:h}];
else if(isLens(o)&&o.power!==undefined){let g=lensGeometry(o),r=g.radius,negative=o.power<0;const half=y=>!Number.isFinite(r)?4:negative?4+y*y/(r+Math.sqrt(r*r-y*y)):4+(h*h-y*y)/(Math.sqrt(r*r-y*y)+Math.sqrt(r*r-h*h));for(let i=0;i<=240;i++){let y=-h+2*h*i/240;p.push({x:half(y),y});}for(let i=240;i>=0;i--){let y=-h+2*h*i/240;p.push({x:-half(y),y});}}
else if(isLens(o)){return polygon(normalizeLens(o));}
else p=[{x:-w,y:-h},{x:w,y:-h},{x:w,y:h},{x:-w,y:h}];let a=rad(o.angle),c=Math.cos(a),s=Math.sin(a);return p.map(v=>({x:o.x+v.x*c-v.y*s,y:o.y+v.x*s+v.y*c}));}
// Edge intersections plus containment work for convex AND concave lenses.
function overlap(a,b){
 const bounds=p=>({x0:Math.min(...p.map(v=>v.x)),x1:Math.max(...p.map(v=>v.x)),y0:Math.min(...p.map(v=>v.y)),y1:Math.max(...p.map(v=>v.y))});let A=bounds(a),B=bounds(b),gap=.25;
 if(A.x1+gap<B.x0||B.x1+gap<A.x0||A.y1+gap<B.y0||B.y1+gap<A.y0)return false;
 const distance=(p,a,b)=>{let v=sub(b,a),l=dot(v,v),t=l?Math.max(0,Math.min(1,dot(sub(p,a),v)/l)):0;return Math.hypot(p.x-a.x-t*v.x,p.y-a.y-t*v.y);};
 for(let i=0;i<a.length;i++)for(let j=0;j<b.length;j++){let p=a[i],q=a[(i+1)%a.length],r=b[j],t=b[(j+1)%b.length];if(Math.max(p.x,q.x)+gap<Math.min(r.x,t.x)||Math.max(r.x,t.x)+gap<Math.min(p.x,q.x)||Math.max(p.y,q.y)+gap<Math.min(r.y,t.y)||Math.max(r.y,t.y)+gap<Math.min(p.y,q.y))continue;
 let v=sub(q,p),w=sub(t,r),den=cross(v,w);if(Math.abs(den)>1e-10){let u=cross(sub(r,p),w)/den,z=cross(sub(r,p),v)/den;if(u>=0&&u<=1&&z>=0&&z<=1)return true;}
 if(Math.min(distance(p,r,t),distance(q,r,t),distance(r,p,q),distance(t,p,q))<=gap)return true;
 }return inside(a[0],b)||inside(b[0],a);
}
function inside(p,poly){let yes=false;for(let i=0,j=poly.length-1;i<poly.length;j=i++){let a=poly[i],b=poly[j];if((a.y>p.y)!=(b.y>p.y)&&p.x<(b.x-a.x)*(p.y-a.y)/(b.y-a.y)+a.x)yes=!yes;}return yes;}
function bend(d,n,n1,n2){if(dot(d,n)>0)n={x:-n.x,y:-n.y};let c=-dot(d,n),r=n1/n2,k=1-r*r*(1-c*c);return k<0?{d:{x:d.x+2*c*n.x,y:d.y+2*c*n.y},tir:true}:{d:{x:r*d.x+(r*c-Math.sqrt(k))*n.x,y:r*d.y+(r*c-Math.sqrt(k))*n.y},tir:false};}
function trace(origin,d,objects,max=40){let p={...origin},segments=[],hits=[],limited=false;const shapes=objects.filter(o=>!['laser','beam'].includes(o.type)).map(o=>({o,poly:polygon(o)}));let medium=shapes.find(s=>s.o.material==='glass'&&inside(p,s.poly));for(let bounce=0;bounce<max;bounce++){let best=null;for(const s of shapes)for(let i=0;i<s.poly.length;i++){let a=s.poly[i],e=sub(s.poly[(i+1)%s.poly.length],a),den=cross(d,e);if(Math.abs(den)<1e-10)continue;let ap=sub(a,p),t=cross(ap,e)/den,u=cross(ap,d)/den;if(t>1e-5&&u>=-1e-8&&u<=1+1e-8&&(!best||t<best.t)){let l=Math.hypot(e.x,e.y);if(l)best={t,s,n:{x:-e.y/l,y:e.x/l}};}}if(!best){segments.push([p,{x:p.x+d.x*6000,y:p.y+d.y*6000}]);break;}let q={x:p.x+d.x*best.t,y:p.y+d.y*best.t};segments.push([p,q]);hits.push({p:q,n:best.n});let o=best.s.o;if(o.material==='absorb')break;if(o.material==='mirror'){let z=dot(d,best.n);d={x:d.x-2*z*best.n.x,y:d.y-2*z*best.n.y};}else{let leaving=medium?.o.id===o.id,r=bend(d,best.n,leaving?o.n:1,leaving?1:o.n);d=r.d;if(!r.tir)medium=leaving?null:best.s;}p={x:q.x+d.x*1e-4,y:q.y+d.y*1e-4};if(bounce===max-1)limited=true;}return {segments,hits,limited};}
// Fresnel power coefficients for lossless dielectrics; s/p components remain separate.
function fresnel(d,n,n1,n2){let c=Math.min(1,Math.abs(dot(d,n))),k=1-(n1/n2)**2*(1-c*c);if(k<0)return {Rs:1,Rp:1,tir:true};if(n1===n2)return {Rs:0,Rp:0,tir:false};let t=Math.sqrt(k),rs=(n1*c-n2*t)/(n1*c+n2*t),rp=(n2*c-n1*t)/(n2*c+n1*t);return {Rs:rs*rs,Rp:rp*rp,tir:false};}
function prepareScene(objects){return objects.filter(o=>!['laser','beam'].includes(o.type)).map(o=>({o,poly:polygon(o)}));}
function traceSplit(origin,d,objects,options={}){
 const shapes=options.shapes||prepareScene(objects),cutoff=options.cutoff??.001,max=options.maxDepth??40,budget=options.budget??256;
 let segments=[],hits=[],escaped=0,absorbed=0,discarded=0,limited=false,branches=0;
 let stack=[{p:{...origin},d,medium:shapes.find(s=>s.o.material==='glass'&&inside(origin,s.poly)),s:.5,powerP:.5,depth:0}];
 const enqueue=(q,d,medium,s,powerP,depth)=>{let power=s+powerP;if(power<cutoff){discarded+=power;return;}stack.push({p:{x:q.x+d.x*1e-4,y:q.y+d.y*1e-4},d,medium,s,powerP,depth});};
 while(stack.length){let ray=stack.pop(),{p,d,medium}=ray,power=ray.s+ray.powerP;if(segments.length>=budget||ray.depth>=max){discarded+=power;limited=true;continue;}
 let best=null;for(const shape of shapes)for(let i=0;i<shape.poly.length;i++){let a=shape.poly[i],e=sub(shape.poly[(i+1)%shape.poly.length],a),den=cross(d,e);if(Math.abs(den)<1e-10)continue;let ap=sub(a,p),t=cross(ap,e)/den,u=cross(ap,d)/den;if(t>1e-5&&u>=-1e-8&&u<=1+1e-8&&(!best||t<best.t)){let l=Math.hypot(e.x,e.y);if(l)best={t,shape,n:{x:-e.y/l,y:e.x/l}};}}
 if(!best){segments.push([p,{x:p.x+d.x*6000,y:p.y+d.y*6000},power]);escaped+=power;continue;}
 let q={x:p.x+d.x*best.t,y:p.y+d.y*best.t},o=best.shape.o,z=dot(d,best.n),reflection={x:d.x-2*z*best.n.x,y:d.y-2*z*best.n.y};segments.push([p,q,power]);hits.push({p:q,n:best.n,power});
 if(o.material==='absorb'){absorbed+=power;continue;}if(o.material==='mirror'){enqueue(q,reflection,medium,ray.s,ray.powerP,ray.depth+1);continue;}
 let leaving=medium?.o.id===o.id,n1=leaving?o.n:1,n2=leaving?1:o.n,f=fresnel(d,best.n,n1,n2);
 if(f.tir){enqueue(q,reflection,medium,ray.s,ray.powerP,ray.depth+1);continue;}
 if(f.Rs+f.Rp>0){branches++;enqueue(q,reflection,medium,ray.s*f.Rs,ray.powerP*f.Rp,ray.depth+1);}
 let transmitted=bend(d,best.n,n1,n2).d;enqueue(q,transmitted,leaving?null:best.shape,ray.s*(1-f.Rs),ray.powerP*(1-f.Rp),ray.depth+1);
 }return {segments,hits,limited,branches,escaped,absorbed,discarded};
}
function snapPoint(p,step){return {x:Math.round(p.x/step)*step,y:Math.round(p.y/step)*step};}
function snappedBox(start,end,step){let a=snapPoint(start,step),b=snapPoint(end,step);if(b.x===a.x)b.x+=step;if(b.y===a.y)b.y+=step;return {x:(a.x+b.x)/2,y:(a.y+b.y)/2,w:Math.abs(b.x-a.x),h:Math.abs(b.y-a.y)};}
root.Optics={fresnel,traceSplit,prepareScene,snapPoint,snappedBox,polygon,overlap,inside,bend,trace,rad,isLens,lensGeometry,lensPower,normalizeLens};if(typeof module!=='undefined')module.exports=root.Optics;
})(typeof window!=='undefined'?window:globalThis);
