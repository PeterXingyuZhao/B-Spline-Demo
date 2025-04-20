import { Mesh, Texture } from '../../../build/three.module.js';

export class TextureHelper extends Mesh {
    texture: Texture;
    type: "TextureHelper";

    constructor(texture: Texture, width?: number, height?: number, depth?: number);

    dispose(): void;
}
