<!--
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
-->
<!DOCTYPE html>
<html>
  <head>
    <title>Codebase Alpha Space Travel</title>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
    <script src="webgl-driver.js" type="text/javascript"></script>
    <link href='https://fonts.googleapis.com/css?family=Josefin Slab' rel='stylesheet'>
    <style>
        body {
            background-color: black;
            color: white;
            font-family: 'Josefin Slab';font-size: 32px;
        }
        canvas.gl {
          position:fixed;
          z-index:-1;
          left:0;
          top:0;
          width:100%;
          height:100%;
        }
    </style>
    <script id="shader-fs" type="x-shader/x-fragment">
      // -----------------------------------------------------------------------
      // BEGIN - Common prelude
      // -----------------------------------------------------------------------
      precision mediump float;

      uniform vec2 iResolution;
      uniform float iTime;
      uniform sampler2D iChannel0;

      varying highp vec2 vTextureCoord;

      void mainImage(out vec4 fragColor, in vec2 fragCoord);

      void main(void) {
        mainImage(gl_FragColor, vTextureCoord*iResolution);
      }
      // -----------------------------------------------------------------------
      // END - Common prelude
      // -----------------------------------------------------------------------

      #define PI              3.141592654
      #define TAU             (2.0*PI)
      #define TIME            iTime
      #define RESOLUTION      iResolution
      #define PERIOD          30.0
      #define FADE            0.0
      
      #define SCA(a)          vec2(sin(a), cos(a))
      #define LESS(a,b,c)     mix(a,b,step(0.,c))
      #define SABS(x,k)       LESS((.5/k)*x*x+k*.5,abs(x),abs(x)-k)
      #define MROT(a)         mat2(cos(a), sin(a), -sin(a), cos(a))
      
      // Set effects below 0, 1 or 2
      #define EFFECT_STRAIGHT 0
      #define EFFECT_CASUAL   0
      #define EFFECT_WARP     0
      
      const mat2 rotSome      = MROT(1.0);
      const vec3 stdGamma     = vec3(2.2);

      float tanh(float x)
    {
      float ex = exp(x);
      float nex = exp(-x);
      float sum = ex + nex;
      return ex/sum - nex/sum;
    }
      
      vec3 hsv2rgb(vec3 c) {
        const vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
        vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
        return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
      }
      
      float hash(in float co) {
        return fract(sin(co*12.9898) * 13758.5453);
      }
      
      float hash(in vec2 co) {
        return fract(sin(dot(co, vec2(12.9898,58.233))) * 13758.5453);
      }
      
      float psin(float a) {
        return 0.5 + 0.5*sin(a);
      }
      
      float pcos(float a) {
        return 0.5 + 0.5*cos(a);
      }
      
      vec2 mod2_1(inout vec2 p) {
        vec2 c = floor(p + 0.5);
        p = fract(p + 0.5) - 0.5;
        return c;
      }
      
      float vnoise(vec2 x) {
        vec2 i = floor(x);
        vec2 w = fract(x);
          
      #if 1
        // quintic interpolation
        vec2 u = w*w*w*(w*(w*6.0-15.0)+10.0);
      #else
        // cubic interpolation
        vec2 u = w*w*(3.0-2.0*w);
      #endif    
      
        float a = hash(i+vec2(0.0,0.0));
        float b = hash(i+vec2(1.0,0.0));
        float c = hash(i+vec2(0.0,1.0));
        float d = hash(i+vec2(1.0,1.0));
          
        float k0 =   a;
        float k1 =   b - a;
        float k2 =   c - a;
        float k3 =   d - c + a - b;
      
        float aa = mix(a, b, u.x);
        float bb = mix(c, d, u.x);
        float cc = mix(aa, bb, u.y);
        
        return k0 + k1*u.x + k2*u.y + k3*u.x*u.y;
      }
      
      vec4 alphaBlendGamma(vec4 back, vec4 front, vec3 gamma) {
        vec3 colb = max(back.xyz, 0.0);
        vec3 colf = max(front.xyz, 0.0);
        colb = pow(colb, gamma);
        colf = pow(colf, gamma);
        vec3 xyz = mix(colb*back.w, colf.xyz, front.w);
        float w = mix(back.w, 1.0, front.w);
        return max(vec4(pow(xyz, 1.0/gamma), w), 0.0);
      }
      
      vec3 alphaBlendGamma(vec3 back, vec4 front, vec3 gamma) {
        vec3 colb = max(back.xyz, 0.0);
        vec3 colf = max(front.xyz, 0.0);;
        
        colb = pow(colb, gamma);
        colf = pow(colf, gamma);
        vec3 xyz = mix(colb, colf.xyz, front.w);
        return pow(xyz, 1.0/gamma);
      }
      
      vec3 offset_0(float z) {
        float a = z;
        vec2 p = vec2(0.0);
        return vec3(p, z);
      }
      
      vec3 offset_1(float z) {
        const float off = 5.0;
        float a = 0.125*z;
        vec2 p = -1.0*(vec2(cos(a), sin(a*sqrt(2.0))) + vec2(cos(a*sqrt(0.75)), sin(a*sqrt(0.5))));
        return vec3(p, z);
      }
      
      vec3 offset(int effect, float z) {
        if (effect == EFFECT_STRAIGHT)    return offset_0(z);
        else if (effect == EFFECT_CASUAL) return offset_1(z);
        else if (effect == EFFECT_WARP)   return offset_0 (z);
        else                              return offset_0 (z);
      }
      
      // Approximate derivate
      vec3 doffset(int effect, float z) {
        float eps = 0.1;
        return 0.5*(offset(effect, z + eps) - offset(effect, z - eps))/eps;
      }
      
      // Approximate 2nd derivate
      vec3 ddoffset(int effect, float z) {
        float eps = 0.1;
        return 0.125*(doffset(effect, z + eps) - doffset(effect, z - eps))/eps;
      }
      
      float globalCloudDensity(vec2 p, float n) {
        p*=0.1;
        float gcd = vnoise(p+10.0*hash(n)+100.0);
        return gcd;
      }
      
      float localCloudDensity(vec2 p, float n) {
        p*=1.0;
        const float aa = -0.45;
        const mat2 pp = 2.03*rotSome;
        float a = 0.5;
        float s = 0.0;
        p += 10.0*hash(n)+100.0;
      
        s += a*vnoise(p); a *= aa; p *= pp;
        s += a*vnoise(p); a *= aa; p *= pp;
        s += a*vnoise(p); a *= aa; p *= pp;
        s += a*vnoise(p); a *= aa; p *= pp;
        s += a*vnoise(p); a *= aa; p *= pp;
        return s*2.1-0.0;
      }
      
      vec4 plane(vec3 ro, vec3 rd, vec3 pp, float aa, float n) {
        vec2 p = pp.xy;
        float z = pp.z;
        float nz = pp.z-ro.z;
      
        float ds = 1E6;
        float r = 0.0;
      
        float gcd = globalCloudDensity(p, n);
      
        float s = mix(1.0, 2.0, gcd);
        vec2 ps = p;
        ps *= s;
        for (int i = 0; i < 5; ++i) {
          ps *= rotSome;
          vec2 ips = ps;
          vec2 ipn = mod2_1(ips);
          float ir = hash(ipn+n*100.0+float(i)*10.0);
          ips = ips -0.3*vec2(ir, fract(-23.0*ir));
          float ids = length(ips)-0.0025;
          if (ids < ds) {
            r = ir;
            ds = min(ds, ids);
          }
        }
      
        float hues = mix(0.6, 0.8, r*r);
        float sats = mix(0.6, 0.0, sqrt(r));
        float bris = mix(0.5, 1.0, r);
        float ts = pow(max(1.0-ds, 0.0), mix(200.0, 100.0, gcd)/sqrt(s));
        vec3 cols = 3.0*hsv2rgb(vec3(hues, sats, bris));
        vec4 cs = vec4(cols, ts);
      
        float cd = gcd*localCloudDensity(p, n);
        float cdo = gcd*localCloudDensity(p+vec2(0.125, 0.25), n);
        const float level0 = 0.0;
        const float level1 = 0.05;
        // Some serious fake shadow & lighting of clouds
        float cli = mix(-0.1, 1.0, 0.5 + 0.5*tanh(10.0*(cd-cdo)));
      
        float huec = (mix(-0.2, 0.05, (cd))+0.05)-0.1*pcos(2.0*pp.z);
        float tc = clamp(cd, 0.0, 1.0);
        float satc = 0.5;
        float bric = 1.0;
        vec3 colc = hsv2rgb(vec3(huec, satc, bric))+cli*0.65;
        vec4 cc = vec4(colc, tc);
        
        vec4 ct = alphaBlendGamma(cs, cc, stdGamma);
        return ct;
      }
      
      vec3 skyColor(vec3 ro, vec3 rd) {
        const vec3 l = normalize(vec3(0.0, 0.0, 1));
        const vec3 baseCol = vec3(0.75, 0.75, 1.0);
        float diff = pow(max(dot(l, rd), 0.0), 20.0);  
        return 1.5*baseCol*pow(max(dot(l, rd), 0.0), 20.0);
      }
      
      vec3 color(vec3 ww, vec3 uu, vec3 vv, vec3 ro, vec2 p, int effect, float effectTime) {
        float lp = length(p);
      
        vec3 rd;
        float planeDist;
        
        if (effect == EFFECT_WARP) {
            float warpf = smoothstep(FADE*2.0, PERIOD, effectTime);
            planeDist = mix(6.0, 0.2, warpf);
            rd = normalize(p.x*uu + p.y*vv + (2.0+mix(0.0, 15.0, warpf)*tanh(lp))*ww);
        } else {
            rd = normalize(p.x*uu + p.y*vv + (2.0+0.75*tanh(lp))*ww);
            planeDist = 3.0;
        }
        
        const int furthest = 6;
        const int fadeFrom = 4; //max((furthest - 2), 0);
      
        float nz = floor(ro.z / planeDist);
      
        vec3 skyCol = skyColor(ro, rd);  
        
        vec3 col = skyCol;
      
        for (int i = furthest; i >= 1 ; --i) {
          float pz = planeDist*nz + planeDist*float(i);
          
          float pd = (pz - ro.z)/rd.z;
          
          if (pd > 0.0) {
            vec3 pp = ro + rd*pd;
         
            float aa = length(pp); //length(dFdy(pp));
      
            vec4 pcol = plane(ro, rd, pp, aa, nz+float(i));
            float nz = pp.z-ro.z;
            float fadeIn = (1.0-smoothstep(planeDist*float(fadeFrom), planeDist*float(furthest), nz));
            float fadeOut = smoothstep(0.0, planeDist*0.25, nz);
            const float rs = 200.0;
            fadeIn *= rs/(rs+pd*pd);
            pcol.xyz = mix(skyCol, pcol.xyz, fadeIn);
            pcol.w *= fadeOut;
      
            col = alphaBlendGamma(col, pcol, stdGamma);
          } else {
            break;
          }
          
        }
        
        return col;
      }
      
      vec3 postProcess(vec3 col, vec2 q)  {
        col=pow(clamp(col,0.0,1.0),vec3(0.75)); 
        col=col*0.6+0.4*col*col*(3.0-2.0*col);
        col=mix(col, vec3(dot(col, vec3(0.33))), -0.4);
        col*=0.5+0.5*pow(19.0*q.x*q.y*(1.0-q.x)*(1.0-q.y),0.7);
        return col;
      }
      
      vec3 effect(vec2 p, vec2 q) {
        float tm = TIME;
        
        float effectTime   = mod(TIME, PERIOD);
        float effectPeriod = floor(TIME/PERIOD);
        
        int effect = int(mod(effectPeriod, 3.0));
        if (effect == EFFECT_WARP) {
          tm = effectTime;
        }
        vec3 ro    = offset(effect, tm);
        vec3 dro   = doffset(effect, tm);
        vec3 ddro  = ddoffset(effect, tm);
      
        vec3 ww = normalize(dro);
        vec3 uu = normalize(cross(normalize(vec3(0.0,1.0,0.0)+ddro), ww));
        vec3 vv = normalize(cross(ww, uu));
        
        vec3 col = color(ww, uu, vv, ro, p, effect, effectTime);
        col = postProcess(col, q);
        
        col *= smoothstep(0.0, FADE*FADE, effectTime);
        col *= smoothstep(0.0, FADE*FADE, PERIOD-effectTime);
        
        return col;
      }
      
      void mainImage(out vec4 fragColor, vec2 fragCoord) {
        vec2 q = fragCoord/RESOLUTION.xy;
        vec2 p = -1. + 2. * q;
        p.x *= RESOLUTION.x/RESOLUTION.y;
        
        vec3 col = effect(p, q);
        
        fragColor = vec4(col, 1.0);
      }
    </script>

    <script id="shader-vs" type="x-shader/x-vertex">
      attribute highp vec3 aVertexPosition;
      attribute highp vec3 aVertexNormal;
      attribute highp vec2 aTextureCoord;

      varying highp vec2 vTextureCoord;
      varying highp vec3 vNormal;

      void main(void) {
        gl_Position   = vec4(aVertexPosition, 1.0);
        vNormal = aVertexNormal;
        vTextureCoord = aTextureCoord;
      }
    </script>
  </head>

  <body onload="start()">
    <canvas id="glcanvas" class="gl">
      Your browser doesn't appear to support the HTML5 <code>&lt;canvas&gt;</code> element.
    </canvas>
  </body>
</html>