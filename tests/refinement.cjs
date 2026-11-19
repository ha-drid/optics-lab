const A=require('node:assert/strict'),O=require('../physics.js');const near=(a,b,e=1e-9)=>A(Math.abs(a-b)<e,`${a} != ${b}`);
let f=O.fresnel({x:1,y:0},{x:-1,y:0},1,1.5);near(f.Rs,.04);near(f.Rp,.04);
let a=Math.atan(1.5);f=O.fresnel({x:Math.cos(a),y:Math.sin(a)},{x:-1,y:0},1,1.5);near(f.Rp,0);A(f.Rs>0);
f=O.fresnel({x:.5,y:Math.sqrt(3)/2},{x:-1,y:0},1.5,1);A(f.tir);near(f.Rs,1);
const box={id:'b',type:'rect',x:0,y:0,w:40,h:400,angle:0,n:1.5,material:'glass'};
let r=O.traceSplit({x:-100,y:0},{x:1,y:0},[box]);A(r.branches>=2);A(r.segments.some(s=>Math.abs(s[2]-.04)<1e-9));A(r.segments.some(s=>Math.abs(s[2]-.9216)<1e-9));near(r.escaped+r.absorbed+r.discarded,1);A(r.discarded<.001);
r=O.traceSplit({x:-100,y:0},{x:1,y:0},[{...box,n:1}]);A.equal(r.branches,0);near(r.escaped,1);
r=O.traceSplit({x:-100,y:0},{x:1,y:0},[{...box,material:'absorb'}]);near(r.absorbed,1);
r=O.traceSplit({x:-100,y:0},{x:1,y:0},[box],{budget:2});A(r.limited);near(r.escaped+r.discarded+r.absorbed,1);A(r.segments.length<=2);
const lens={...box,type:'diverging',h:240,power:3};let l=O.normalizeLens(lens);A.equal(l.type,'diverging');A.equal(l.power,-3);l=O.normalizeLens({...l,power:0});A.equal(l.type,'diverging');l=O.normalizeLens({...l,power:2});A.equal(l.power,-2);A.equal(O.normalizeLens({...l,type:'lens'}).power,2);
for(let step of [10,25,50,100]){let b=O.snappedBox({x:14,y:13},{x:115,y:94},step);A.equal(b.w%step,0);A.equal(b.h%step,0);near((b.x-b.w/2)%step,0);near((b.y-b.h/2)%step,0);}
console.log('Passed: Fresnel 4%/96%, Brewster, TIR, branching, energy conservation, cutoff/budget, fixed lens type, snapped block corners');
