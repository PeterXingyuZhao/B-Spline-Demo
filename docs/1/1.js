/*jshint esversion: 6 */
// @ts-check

import * as T from "../../libs/CS559-Three/build/three.module.js";
import { GrWorld } from "../../libs/CS559-Framework/GrWorld.js";
import { OrbitControls } from "../../libs/CS559-Three/examples/jsm/controls/OrbitControls.js";



let div = document.getElementById("div1");
let world = new GrWorld({ groundplanesize: 20, where: div, renderparams: {preserveDrawingBuffer:true}, id:"canvas" });

world.go();


/* non-framework version 

let renderer = new T.WebGLRenderer({preserveDrawingBuffer:true});
renderer.setSize(500, 500);
document.getElementById("div1").appendChild(renderer.domElement);
renderer.domElement.id = "canvas";
*/ 

