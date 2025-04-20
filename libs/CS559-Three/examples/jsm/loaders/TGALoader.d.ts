import { DataTexture, DataTextureLoader, LoadingManager } from '../../../build/three.module.js';

export class TGALoader extends DataTextureLoader {
    constructor(manager?: LoadingManager);

    parse(data: ArrayBuffer): DataTexture;
}
