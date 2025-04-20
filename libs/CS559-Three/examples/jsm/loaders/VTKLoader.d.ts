import { BufferGeometry, Loader, LoadingManager } from '../../../build/three.module.js';

export class VTKLoader extends Loader<BufferGeometry> {
    constructor(manager?: LoadingManager);

    parse(data: ArrayBuffer | string, path: string): BufferGeometry;
}
