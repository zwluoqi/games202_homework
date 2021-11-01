class HW2Material extends Material {

    constructor( vertexShader, fragmentShader) {


        super({
            //'uLightMVP': { type: 'matrix4fv', value: lightMVP }
            
        }, ['aPrecomputeLT'], vertexShader, fragmentShader, null);
    }
}

async function buildPRTMaterial(vertexPath, fragmentPath) {
    

    let vertexShader = await getShaderString(vertexPath);
    let fragmentShader = await getShaderString(fragmentPath);

    return new HW2Material(vertexShader, fragmentShader);

}