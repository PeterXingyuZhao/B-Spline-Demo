/*jshint esversion: 11 */
// @ts-check

/**
 * CS559 3D World Framework Code
 *
 * Simple, automatic UI from an object with properly declared parameters
 * 
 * @module AutoUI 
 * */

// we need to have the BaseClass definition
import { GrObject } from "./GrObject.js";
// we need to import the module to get its typedefs for the type checker
import * as InputHelpers from "../CS559/inputHelpers.js";
import { GUI } from '../CS559-Three/examples/jsm/libs/lil-gui.module.min.js';

/**
 * This is a "global" variable - if panels are placed without a where,
 * we make a DIV (the "panel panel") and put them in there - this way
 * we can get floating
 * 
 * @type{HTMLElement}
 */
let panelPanel;
// since exports are read only, always access it by a function that will make it
export function panel() {
  if (!panelPanel) {
    panelPanel = InputHelpers.makeFlexDiv();
    panelPanel.id = "panel-panel";
  }
  return panelPanel;
}

export class AutoUI {
  /**
   * Create a UI panel for a GrObject
   * goes through the parameters and makes a slider for each
   * also defines a callback for those sliders that calls the
   * object's update function.
   *
   * This does place the panel into the DOM (onto the web page)
   * using `insertElement` in the CS559 helper library. The place
   * it is placed is controlled the `where` parameter. By default,
   * it goes at the end of the DOM. However, you can pass it a DOM
   * element to be placed inside (or some other choices as well).
   *
   * @param {GrObject} object
   * @param {number} [width=300]
   * @param {InputHelpers.WhereSpec} [where] - where to place the panel in the DOM (at the end of the page by default)
   * @param {boolean} adjusted - whether adjust the slider length according to the label length
   * @param {string} display - align type of the label and slider
   */
  constructor(object, width = 300, where = undefined, widthdiv = 1, adjusted = false, display = "") {
    const self = this;
    this.object = object;

    // Create the GUI using lil-gui
    let element = document.getElementById("gui");
    let gui;
    if (!element) {
      if (adjusted) gui = new GUI({ title: "AutoUI" });
      else gui = new GUI({ width: width / widthdiv, title: "AutoUI" });
      gui.domElement.id = "gui";
      gui.domElement.gui = gui;
    }
    else gui = element.gui;
    const folder = gui.addFolder(object.name);
    object.params.forEach(function (param) {
      if (object.values) object.values[param.name] = param.initial;
      else object.values = { [param.name]: param.initial };
      folder.add(object.values, param.name, param.min, param.max, param.step || Math.max((param.max - param.min) / 30, 1)).onChange(function () {
        object.update(folder.controllers.map(c => c.getValue()));
      });
    });
    folder.close();

    if (display) {
      /* if no where is provided, put it at the end of the panel panel - assuming there is one */
      if (!where) {
        where = panel();
      }

      // Create the sliders using the CS559 inputHelpers
      this.div = InputHelpers.makeBoxDiv({ width: width, flex: widthdiv > 1 }, where);
      InputHelpers.makeHead(object.name, this.div, { tight: true });
      if (widthdiv > 1) InputHelpers.makeFlexBreak(this.div);

      this.sliders = object.params.map(function (param) {
        const slider = new InputHelpers.LabelSlider(param.name, {
          where: self.div,
          width: (width / widthdiv) - 20,
          min: param.min,
          max: param.max,
          step: param.step ?? ((param.max - param.min) / 30),
          initial: param.initial,
          id: object.name + "-" + param.name,
          adjusted: adjusted,
          display: display,
        });
        return slider;
      });

      this.sliders.forEach(function (sl) {
        sl.oninput = function () {
          self.update();
        };
      });

      this.update();
    }
  }
  update() {
    const vals = this.sliders.map(sl => Number(sl.value()));
    this.object.update(vals);
  }

  /**
   *
   * @param {number | string} param
   * @param {number} value
   */
  set(param, value) {
    if (typeof param === "string") {
      for (let i = 0; i < this.object.params.length; i++) {
        if (param == this.object.params[i].name) {
          this.sliders[i].set(Number(value));
          return;
        }
      }
      throw `Bad parameter ${param} to set`;
    } else {
      this.sliders[param].set(Number(value));
    }
  }
}