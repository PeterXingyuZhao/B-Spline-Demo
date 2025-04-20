import { Loader, Scene } from '../../../build/three.module.js';

export interface Collada {
    kinematics: object;
    library: object;
    scene: Scene;
}

export class ColladaLoader extends Loader<Collada> {
    parse(text: string, path: string): Collada;
}
