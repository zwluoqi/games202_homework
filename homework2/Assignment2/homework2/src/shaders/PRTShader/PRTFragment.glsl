#ifdef GL_ES
precision mediump float;
#endif
//uniform sampler2D uSampler;

//varying highp vec2 vTextureCoord;
varying highp vec3 vNormal;
varying highp vec3 vColor;


void main(void){

      //vec3 color = texture2D(uSampler, vTextureCoord).rgb;
      gl_FragColor = vec4(vColor, 1.0);

      
    //  gl_FragColor = vec4(vColor.g,0.0,0.0, 1.0);
    //gl_FragColor = vec4(1.0,1.0,1.0, 1.0);
}