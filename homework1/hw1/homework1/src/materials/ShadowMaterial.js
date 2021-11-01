class ShadowMaterial extends Material {

    constructor(light, vertexShader, fragmentShader) {
        // let lightMVP = light.CalcLightMVP(translate, scale);

        super({
            // 'uLightMVP': { type: 'matrix4fv', value: lightMVP }
            // 'uLightPos':{},
            // 'uLightMVP':{},
        }, [], vertexShader, fragmentShader, light.fbo);
    }
}

async function buildShadowMaterial(light, vertexPath, fragmentPath) {


    let vertexShader = await getShaderString(vertexPath);
    let fragmentShader = await getShaderString(fragmentPath);

    return new ShadowMaterial(light, vertexShader, fragmentShader);

}