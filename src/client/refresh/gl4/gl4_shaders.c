/*
 * Copyright (C) 1997-2001 Id Software, Inc.
 * Copyright (C) 2016-2017 Daniel Gibson
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 *
 * =======================================================================
 *
 * OpenGL4 refresher: Handling shaders
 *
 *   *** PS1-style vertex rounding ***
 *
 *  - 3D vertex shaders snap their vertex positions to a low-resolution
 *    grid before the matrix multiplication, giving the classic
 *    "wobbly" PlayStation 1 look.
 *  - No integer attributes, no VAO changes, no glVertexAttribIPointer.
 *  - A float uniform 'cameraPos' can optionally be subtracted before
 *    rounding to keep precision on large maps.
 *  - The grid size is exposed as a float uniform 'ps1GridSize' so the
 *    engine can change it at runtime (smaller = chunkier).
 *
 * =======================================================================
 */

#include "header/local.h"

// TODO: remove eprintf() usage
#define eprintf(...)  R_Printf(PRINT_ALL, __VA_ARGS__)

// Set to 1 to snap AFTER projection (screen space, most PS1-accurate).
// Set to 0 to snap world-space vertex positions (usually looks better
// in a modern engine and still gives the retro wobble).
#define PS1_SCREEN_SPACE_SNAP 1


static GLuint
CompileShader(GLenum shaderType, const char* shaderSrc, const char* shaderSrc2)
{
	GLuint shader = glCreateShader(shaderType);
	const char* version = glshader_version(gl4config.major_version, gl4config.minor_version);
	const char* sources[3] = { version, shaderSrc, shaderSrc2 };
	int numSources = shaderSrc2 != NULL ? 3 : 2;

	glShaderSource(shader, numSources, sources, NULL);
	glCompileShader(shader);
	GLint status;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
	if (status != GL_TRUE)
	{
		char buf[2048];
		char* bufPtr = buf;
		int bufLen = sizeof(buf);
		GLint infoLogLength;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogLength);
		if (infoLogLength >= bufLen)
		{
			bufPtr = malloc(infoLogLength+1);
			bufLen = infoLogLength+1;
			if (bufPtr == NULL)
			{
				bufPtr = buf;
				bufLen = sizeof(buf);
				eprintf("WARN: In CompileShader(), malloc(%d) failed!\n", infoLogLength+1);
			}
		}

		glGetShaderInfoLog(shader, bufLen, NULL, bufPtr);

		const char* shaderTypeStr = "";
		switch(shaderType)
		{
			case GL_VERTEX_SHADER:   shaderTypeStr = "vertex"; break;
			case GL_FRAGMENT_SHADER: shaderTypeStr = "fragment"; break;
			case GL_COMPUTE_SHADER:  shaderTypeStr = "compute"; break;
		}
		eprintf("ERROR: Compiling %s Shader failed: %s\n", shaderTypeStr, bufPtr);
		glDeleteShader(shader);

		if (bufPtr != buf)  free(bufPtr);

		return 0;
	}

	return shader;
}

static GLuint
CreateShaderProgram(int numShaders, const GLuint* shaders)
{
	int i=0;
	GLuint shaderProgram = glCreateProgram();
	if (shaderProgram == 0)
	{
		eprintf("ERROR: Couldn't create a new Shader Program!\n");
		return 0;
	}

	for (i=0; i<numShaders; ++i)
	{
		glAttachShader(shaderProgram, shaders[i]);
	}

	// make sure all shaders use the same attribute locations for common attributes
	// (so the same VAO can easily be used with different shaders)
	glBindAttribLocation(shaderProgram, GL4_ATTRIB_POSITION, "position");
	glBindAttribLocation(shaderProgram, GL4_ATTRIB_TEXCOORD, "texCoord");
	glBindAttribLocation(shaderProgram, GL4_ATTRIB_LMTEXCOORD, "lmTexCoord");
	glBindAttribLocation(shaderProgram, GL4_ATTRIB_COLOR, "vertColor");
	glBindAttribLocation(shaderProgram, GL4_ATTRIB_NORMAL, "normal");
	glBindAttribLocation(shaderProgram, GL4_ATTRIB_LIGHTFLAGS, "lightFlags");

	// the following line is not necessary/implicit (as there's only one output)
	// glBindFragDataLocation(shaderProgram, 0, "outColor"); XXX would this even be here?

	glLinkProgram(shaderProgram);

	GLint status;
	glGetProgramiv(shaderProgram, GL_LINK_STATUS, &status);
	if (status != GL_TRUE)
	{
		char buf[2048];
		char* bufPtr = buf;
		int bufLen = sizeof(buf);
		GLint infoLogLength;
		glGetProgramiv(shaderProgram, GL_INFO_LOG_LENGTH, &infoLogLength);
		if (infoLogLength >= bufLen)
		{
			bufPtr = malloc(infoLogLength+1);
			bufLen = infoLogLength+1;
			if (bufPtr == NULL)
			{
				bufPtr = buf;
				bufLen = sizeof(buf);
				eprintf("WARN: In CreateShaderProgram(), malloc(%d) failed!\n", infoLogLength+1);
			}
		}

		glGetProgramInfoLog(shaderProgram, bufLen, NULL, bufPtr);

		eprintf("ERROR: Linking shader program failed: %s\n", bufPtr);

		glDeleteProgram(shaderProgram);

		if (bufPtr != buf)  free(bufPtr);

		return 0;
	}

	for (i = 0; i < numShaders; ++i)
	{
		// after linking, they don't need to be attached anymore.
		// no idea  why they even are, if they don't have to..
		glDetachShader(shaderProgram, shaders[i]);
	}

	return shaderProgram;
}

#define MULTILINE_STRING(...) #__VA_ARGS__

// ############## shaders for 2D rendering (HUD, menus, console, videos, ..) #####################
// NOTE: 2D shaders work in screen space and are NOT affected by PS1 rounding.

static const char* vertexSrc2D = MULTILINE_STRING(

		in vec2 position; // GL4_ATTRIB_POSITION
		in vec2 texCoord; // GL4_ATTRIB_TEXCOORD

		// for UBO shared between 2D shaders
		layout (std140) uniform uni2D
		{
			mat4 trans;
		};

		noperspective out vec2 passTexCoord;

		void main()
		{
			gl_Position = trans * vec4(position, 0.0, 1.0);
			passTexCoord = texCoord;
		}
);

static const char* fragmentSrc2D = MULTILINE_STRING(

		noperspective in vec2 passTexCoord;

		// for UBO shared between all shaders (incl. 2D)
		layout (std140) uniform uniCommon
		{
			float gamma;
			float intensity;
			float intensity2D; // for HUD, menu etc

			vec4 color;
		};

		uniform sampler2D tex;

		out vec4 outColor;

		void main()
		{
			vec4 texel = texture(tex, passTexCoord);
			// the gl1 renderer used glAlphaFunc(GL_GREATER, 0.666);
			// and glEnable(GL_ALPHA_TEST); for 2D rendering
			// this should do the same
			if (texel.a <= 0.666)
				discard;

			// apply gamma correction and intensity
			texel.rgb *= intensity2D;
			outColor.rgb = pow(texel.rgb, vec3(gamma));
			outColor.a = texel.a; // I think alpha shouldn't be modified by gamma and intensity
		}
);

// like fragmentSrc2D, but also multiplies by color uniform for tinting (e.g. crosshair color)
static const char* fragmentSrc2Dtinted = MULTILINE_STRING(

		noperspective in vec2 passTexCoord;

		// for UBO shared between all shaders (incl. 2D)
		layout (std140) uniform uniCommon
		{
			float gamma;
			float intensity;
			float intensity2D; // for HUD, menu etc

			vec4 color;
		};

		uniform sampler2D tex;

		out vec4 outColor;

		void main()
		{
			vec4 texel = texture(tex, passTexCoord);

			if (texel.a <= 0.666)
				discard;

			// apply color tint
			texel.rgb *= color.rgb;

			// apply gamma correction and intensity
			texel.rgb *= intensity2D;
			outColor.rgb = pow(texel.rgb, vec3(gamma));
			outColor.a = texel.a;
		}
);

static const char* fragmentSrc2Dpostprocess = MULTILINE_STRING(
		noperspective in vec2 passTexCoord;

		uniform sampler2D tex;
		uniform vec4 v_blend;

		out vec4 outColor;

		void main()
		{
			// no gamma or intensity here, it has been applied before
			// (this is just for postprocessing)
			vec4 res = texture(tex, passTexCoord);
			// apply the v_blend, usually blended as a colored quad with:
			res.rgb = v_blend.a * v_blend.rgb + (1.0 - v_blend.a)*res.rgb;
			outColor =  res;
		}
);

static const char* fragmentSrc2DpostprocessWater = MULTILINE_STRING(
		noperspective in vec2 passTexCoord;

		const float PI = 3.14159265358979323846;

		uniform sampler2D tex;

		uniform float time;
		uniform vec4 v_blend;

		out vec4 outColor;

		void main()
		{
			vec2 uv = passTexCoord;

			// warping based on ref_vk
			float sx = 1.0 - abs(0.5 - uv.x) * 2.0;
			float sy = 1.0 - abs(0.5 - uv.y) * 2.0;
			float xShift = 2.0 * time + uv.y * PI * 10.0;
			float yShift = 2.0 * time + uv.x * PI * 10.0;
			vec2 distortion = vec2(sin(xShift) * sx, sin(yShift) * sy) * 0.00666;

			uv += distortion;
			uv = clamp(uv, vec2(0.0, 0.0), vec2(1.0, 1.0));

			// no gamma or intensity here, it has been applied before
			// (this is just for postprocessing)
			vec4 res = texture(tex, uv);

			// apply the v_blend, usually blended as a colored quad with:
			res.rgb = v_blend.a * v_blend.rgb + (1.0 - v_blend.a) * res.rgb;
			outColor =  res;
		}
);

// 2D color only rendering, GL4_Draw_Fill(), GL4_Draw_FadeScreen()
static const char* vertexSrc2Dcolor = MULTILINE_STRING(

		in vec2 position; // GL4_ATTRIB_POSITION

		// for UBO shared between 2D shaders
		layout (std140) uniform uni2D
		{
			mat4 trans;
		};

		void main()
		{
			gl_Position = trans * vec4(position, 0.0, 1.0);
		}
);

static const char* fragmentSrc2Dcolor = MULTILINE_STRING(

		// for UBO shared between all shaders (incl. 2D)
		layout (std140) uniform uniCommon
		{
			float gamma;
			float intensity;
			float intensity2D; // for HUD, menus etc

			vec4 color;
		};

		out vec4 outColor;

		void main()
		{
			vec3 col = color.rgb * intensity2D;
			outColor.rgb = pow(col, vec3(gamma));
			outColor.a = color.a;
		}
);

// ############## shaders for 3D rendering #####################
// NOTE: 'position' is a normal float vec3 (VAO unchanged).
//       PS1 rounding is done inside each vertex shader's main().

static const char* vertexCommon3D = MULTILINE_STRING(

		in vec3 position;   // GL4_ATTRIB_POSITION (float, unchanged)
		in vec2 texCoord;   // GL4_ATTRIB_TEXCOORD
		in vec2 lmTexCoord; // GL4_ATTRIB_LMTEXCOORD
		in vec4 vertColor;  // GL4_ATTRIB_COLOR
		in vec3 normal;     // GL4_ATTRIB_NORMAL
		in uint lightFlags; // GL4_ATTRIB_LIGHTFLAGS

		noperspective out vec2 passTexCoord;

		// Camera position in world space, used only when
		// subtracting to keep float precision on large maps.
		// Set from CPU via GL4_SetCameraPosition(); harmless if left at 0.
		uniform vec3 cameraPos;

		// PS1 vertex grid size. Smaller = chunkier wobble.
		// Default is set from CPU via GL4_SetPS1Grid().
		uniform float ps1GridSize;

#if PS1_SCREEN_SPACE_SNAP
		// Screen resolution used for snapping in clip space.
		// Adjust to your internal render resolution.
		const vec2 ps1ScreenRes = vec2(320.0, 240.0);
#endif

		// Snap a world-space position to the PS1 grid.
		vec3 ps1RoundWorld(vec3 p)
		{
			float g = max(ps1GridSize, 1e-6);
			return floor(p / g + 0.5) * g;
		}

		// for UBO shared between all 3D shaders
		layout (std140) uniform uni3D
		{
			mat4 transProjView;
			mat4 transModel;

			/* Fog parameters */
			vec4 fogColor; // RGB + density in .w
			vec4 heightfog_start; // RGB + start distance in .w
			vec4 heightfog_end; // RGB + end distance in .w

			float sscroll; // for SURF_FLOWING
			float tscroll; // for SURF_FLOWING
			float time;
			float alpha;
			float overbrightbits;
			float particleFadeFactor;
			float lightScaleForTurb; // surfaces with SURF_DRAWTURB (water, lava) don't have lightmaps, use this instead

			float heightfog_density;
			float heightfog_falloff;
			// AMDs legacy windows driver needs this, otherwise uni3D has wrong non std140 size, round up to 16 bytes?
			float _std140_pad1;
			float _std140_pad2;
			float _std140_pad3;
		};
);

static const char* fragmentCommon3D = MULTILINE_STRING(

		noperspective in vec2 passTexCoord;

		out vec4 outColor;

		// for UBO shared between all shaders (incl. 2D)
		layout (std140) uniform uniCommon
		{
			float gamma; // this is 1.0/vid_gamma
			float intensity;
			float intensity2D; // for HUD, menus etc

			vec4 color; // really?
		};
		// for UBO shared between all 3D shaders
		layout (std140) uniform uni3D
		{
			mat4 transProjView;
			mat4 transModel;

			/* Fog parameters */
			vec4 fogColor; // RGB + density in .w
			vec4 heightfog_start; // RGB + start distance in .w
			vec4 heightfog_end; // RGB + end distance in .w

			float sscroll; // for SURF_FLOWING
			float tscroll; // for SURF_FLOWING
			float time;
			float alpha;
			float overbrightbits;
			float particleFadeFactor;
			float lightScaleForTurb; // surfaces with SURF_DRAWTURB (water, lava) don't have lightmaps, use this instead

			float heightfog_density;
			float heightfog_falloff;
			// AMDs legacy windows driver needs this, otherwise uni3D has wrong non std140 size, round up to 16 bytes?
			float _std140_pad1;
			float _std140_pad2;
			float _std140_pad3;
		};
);

static const char* vertexSrc3D = MULTILINE_STRING(

		// it gets attributes and uniforms from vertexCommon3D

		void main()
		{
			vec3 relPos = position - cameraPos;
			vec3 roundedPos = ps1RoundWorld(relPos);

			passTexCoord = texCoord;
			gl_Position = transProjView * transModel * vec4(roundedPos, 1.0);

#if PS1_SCREEN_SPACE_SNAP
			gl_Position.xy = floor(gl_Position.xy * ps1ScreenRes * 0.5)
			                 / (ps1ScreenRes * 0.5);
#endif
		}
);

static const char* vertexSrc3Dflow = MULTILINE_STRING(

		// it gets attributes and uniforms from vertexCommon3D

		void main()
		{
			vec3 relPos = position - cameraPos;
			vec3 roundedPos = ps1RoundWorld(relPos);

			passTexCoord = texCoord + vec2(sscroll, tscroll);
			gl_Position = transProjView * transModel * vec4(roundedPos, 1.0);

#if PS1_SCREEN_SPACE_SNAP
			gl_Position.xy = floor(gl_Position.xy * ps1ScreenRes * 0.5)
			                 / (ps1ScreenRes * 0.5);
#endif
		}
);

static const char* vertexSrc3Dlm = MULTILINE_STRING(

		// it gets attributes and uniforms from vertexCommon3D

		noperspective out vec2 passLMcoord;
		noperspective out vec3 passWorldCoord;
		noperspective out vec3 passNormal;
		flat out uint passLightFlags;

		void main()
		{
			vec3 relPos = position - cameraPos;
			vec3 roundedPos = ps1RoundWorld(relPos);

			passTexCoord = texCoord;
			passLMcoord = lmTexCoord;
			vec4 worldCoord = transModel * vec4(roundedPos, 1.0);
			// NOTE: camera-relative world coord. dynLights origins should
			// be camera-relative too if you use GL4_SetCameraPosition().
			passWorldCoord = worldCoord.xyz;
			vec4 worldNormal = transModel * vec4(normal, 0.0f);
			passNormal = normalize(worldNormal.xyz);
			passLightFlags = lightFlags;

			gl_Position = transProjView * worldCoord;

#if PS1_SCREEN_SPACE_SNAP
			gl_Position.xy = floor(gl_Position.xy * ps1ScreenRes * 0.5)
			                 / (ps1ScreenRes * 0.5);
#endif
		}
);

static const char* vertexSrc3DlmFlow = MULTILINE_STRING(

		// it gets attributes and uniforms from vertexCommon3D

		noperspective out vec2 passLMcoord;
		noperspective out vec3 passWorldCoord;
		noperspective out vec3 passNormal;
		flat out uint passLightFlags;

		void main()
		{
			vec3 relPos = position - cameraPos;
			vec3 roundedPos = ps1RoundWorld(relPos);

			passTexCoord = texCoord + vec2(sscroll, tscroll);
			passLMcoord = lmTexCoord;
			vec4 worldCoord = transModel * vec4(roundedPos, 1.0);
			passWorldCoord = worldCoord.xyz;
			vec4 worldNormal = transModel * vec4(normal, 0.0f);
			passNormal = normalize(worldNormal.xyz);
			passLightFlags = lightFlags;

			gl_Position = transProjView * worldCoord;

#if PS1_SCREEN_SPACE_SNAP
			gl_Position.xy = floor(gl_Position.xy * ps1ScreenRes * 0.5)
			                 / (ps1ScreenRes * 0.5);
#endif
		}
);

static const char* fragmentSrc3D = MULTILINE_STRING(

		// it gets attributes and uniforms from fragmentCommon3D

		uniform sampler2D tex;

		void main()
		{
			vec4 texel = texture(tex, passTexCoord);

			// apply intensity and gamma
			texel.rgb *= intensity;
			outColor.rgb = pow(texel.rgb, vec3(gamma));
			outColor.a = texel.a*alpha; // I think alpha shouldn't be modified by gamma and intensity

			// Apply global fog if enabled (density > 0)
			if (fogColor.w > 0.0)
			{
				float depth = gl_FragCoord.z / gl_FragCoord.w;
				float d = fogColor.w * depth;
				float fogFactor = 1.0 - exp(-(d * d)); // quadratic exponential falloff
				outColor.rgb = mix(outColor.rgb, fogColor.rgb, fogFactor);
			}
		}
);

static const char* fragmentSrc3Dwater = MULTILINE_STRING(

		// it gets attributes and uniforms from fragmentCommon3D

		uniform sampler2D tex;

		void main()
		{
			vec2 tc = passTexCoord;
			tc.s += sin( passTexCoord.t*0.125 + time ) * 4.0;
			tc.s += sscroll;
			tc.t += sin( passTexCoord.s*0.125 + time ) * 4.0;
			tc.s += tscroll;
			tc *= 1.0/64.0; // do this last

			vec4 texel = texture(tex, tc);

			// apply intensity and gamma
			texel.rgb *= intensity * lightScaleForTurb;
			outColor.rgb = pow(texel.rgb, vec3(gamma));
			outColor.a = texel.a*alpha; // I think alpha shouldn't be modified by gamma and intensity

			// Apply global fog if enabled (density > 0)
			if (fogColor.w > 0.0)
			{
				float depth = gl_FragCoord.z / gl_FragCoord.w;
				float d = fogColor.w * depth;
				float fogFactor = 1.0 - exp(-(d * d)); // quadratic exponential falloff
				outColor.rgb = mix(outColor.rgb, fogColor.rgb, fogFactor);
			}
		}
);

static const char* fragmentSrc3Dlm = MULTILINE_STRING(

		// it gets attributes and uniforms from fragmentCommon3D

		struct DynLight { // gl4UniDynLight in C
			vec3 lightOrigin; // NOTE: if you use GL4_SetCameraPosition(),
			                  // upload this already camera-relative.
			float _pad;
			//vec3 lightColor;
			//float lightIntensity;
			vec4 lightColor; // .a is intensity; this way it also works on OSX...
			// (otherwise lightIntensity always contained 1 there)
		};

		layout (std140) uniform uniLights
		{
			DynLight dynLights[32];
			uint numDynLights;
			uint _pad1; uint _pad2; uint _pad3; // FFS, AMD!
		};

		uniform sampler2D tex;

		uniform sampler2D lightmap0;
		uniform sampler2D lightmap1;
		uniform sampler2D lightmap2;
		uniform sampler2D lightmap3;

		uniform vec4 lmScales[4];

		noperspective in vec2 passLMcoord;
		noperspective in vec3 passWorldCoord; // camera-relative if cameraPos is set
		noperspective in vec3 passNormal;
		flat in uint passLightFlags;

		void main()
		{
			vec4 texel = texture(tex, passTexCoord);

			// apply intensity
			texel.rgb *= intensity;

			// apply lightmap
			vec4 lmTex = texture(lightmap0, passLMcoord) * lmScales[0];
			lmTex     += texture(lightmap1, passLMcoord) * lmScales[1];
			lmTex     += texture(lightmap2, passLMcoord) * lmScales[2];
			lmTex     += texture(lightmap3, passLMcoord) * lmScales[3];

			if (passLightFlags != 0u)
			{
				// TODO: or is hardcoding 32 better?
				for (uint i=0u; i<numDynLights; ++i)
				{
					// dyn light number i does not affect this plane, just skip it
					if ((passLightFlags & (1u << i)) == 0u)  continue;

					float intens = dynLights[i].lightColor.a;

					vec3 lightToPos = dynLights[i].lightOrigin - passWorldCoord;
					float distLightToPos = length(lightToPos);
					float fact = max(0.0, intens - distLightToPos - 52.0);

					// move the light source a bit further above the surface
					lightToPos += passNormal*32.0;

					// also factor in angle between light and point on surface
					fact *= max(0.0, dot(passNormal, normalize(lightToPos)));

					lmTex.rgb += dynLights[i].lightColor.rgb * fact * (1.0/256.0);
				}
			}

			lmTex.rgb *= overbrightbits;
			outColor = lmTex*texel;
			outColor.rgb = pow(outColor.rgb, vec3(gamma)); // apply gamma correction to result

			// Apply fog if enabled
			if (fogColor.w > 0.0)
			{
				float depth = gl_FragCoord.z / gl_FragCoord.w;
				float d = fogColor.w * depth;
				float fogFactor = 1.0 - exp(-(d * d));
				outColor.rgb = mix(outColor.rgb, fogColor.rgb, fogFactor);
			}

			outColor.a = 1.0; // lightmaps aren't used with translucent surfaces
		}
);

static const char* fragmentSrc3DlmNoColor = MULTILINE_STRING(

		// it gets attributes and uniforms from fragmentCommon3D

		struct DynLight { // gl4UniDynLight in C
			vec3 lightOrigin;
			float _pad;
			vec4 lightColor; // .a is intensity
		};

		layout (std140) uniform uniLights
		{
			DynLight dynLights[32];
			uint numDynLights;
			uint _pad1; uint _pad2; uint _pad3; // FFS, AMD!
		};

		uniform sampler2D tex;

		uniform sampler2D lightmap0;
		uniform sampler2D lightmap1;
		uniform sampler2D lightmap2;
		uniform sampler2D lightmap3;

		uniform vec4 lmScales[4];

		noperspective in vec2 passLMcoord;
		noperspective in vec3 passWorldCoord;
		noperspective in vec3 passNormal;
		flat in uint passLightFlags;

		void main()
		{
			vec4 texel = texture(tex, passTexCoord);

			// apply intensity
			texel.rgb *= intensity;

			// apply lightmap
			vec4 lmTex = texture(lightmap0, passLMcoord) * lmScales[0];
			lmTex     += texture(lightmap1, passLMcoord) * lmScales[1];
			lmTex     += texture(lightmap2, passLMcoord) * lmScales[2];
			lmTex     += texture(lightmap3, passLMcoord) * lmScales[3];

			if (passLightFlags != 0u)
			{
				for (uint i=0u; i<numDynLights; ++i)
				{
					if ((passLightFlags & (1u << i)) == 0u)  continue;

					float intens = dynLights[i].lightColor.a;

					vec3 lightToPos = dynLights[i].lightOrigin - passWorldCoord;
					float distLightToPos = length(lightToPos);
					float fact = max(0.0, intens - distLightToPos - 52.0);

					lightToPos += passNormal*32.0;

					fact *= max(0.0, dot(passNormal, normalize(lightToPos)));

					lmTex.rgb += dynLights[i].lightColor.rgb * fact * (1.0/256.0);
				}
			}

			// turn lightcolor into grey for gl4_colorlight 0
			lmTex.rgb = vec3(0.333 * (lmTex.r+lmTex.g+lmTex.b));

			lmTex.rgb *= overbrightbits;
			outColor = lmTex*texel;
			outColor.rgb = pow(outColor.rgb, vec3(gamma)); // apply gamma correction to result

			// Apply fog if enabled
			if (fogColor.w > 0.0)
			{
				float depth = gl_FragCoord.z / gl_FragCoord.w;
				float d = fogColor.w * depth;
				float fogFactor = 1.0 - exp(-(d * d));
				outColor.rgb = mix(outColor.rgb, fogColor.rgb, fogFactor);
			}

			outColor.a = 1; // lightmaps aren't used with translucent surfaces
		}
);

static const char* fragmentSrc3Dcolor = MULTILINE_STRING(

		// it gets attributes and uniforms from fragmentCommon3D

		void main()
		{
			vec4 texel = color;

			// apply gamma correction and intensity
			outColor.rgb = pow(texel.rgb, vec3(gamma));

			// Apply fog if enabled
			if (fogColor.w > 0.0)
			{
				float depth = gl_FragCoord.z / gl_FragCoord.w;
				float d = fogColor.w * depth;
				float fogFactor = 1.0 - exp(-(d * d));
				outColor.rgb = mix(outColor.rgb, fogColor.rgb, fogFactor);
			}

			outColor.a = texel.a*alpha; // I think alpha shouldn't be modified by gamma and intensity
		}
);

static const char* fragmentSrc3Dsky = MULTILINE_STRING(

		// it gets attributes and uniforms from fragmentCommon3D

		uniform sampler2D tex;

		void main()
		{
			vec4 texel = texture(tex, passTexCoord);

			// apply gamma correction
			outColor.rgb = pow(texel.rgb, vec3(gamma));

			// Apply fog if enabled
			if (fogColor.w > 0.0)
			{
				float depth = gl_FragCoord.z / gl_FragCoord.w;
				float d = fogColor.w * depth;
				float fogFactor = 1.0 - exp(-(d * d));
				outColor.rgb = mix(outColor.rgb, fogColor.rgb, fogFactor);
			}

			outColor.a = texel.a*alpha; // I think alpha shouldn't be modified by gamma and intensity
		}
);

static const char* fragmentSrc3Dsprite = MULTILINE_STRING(

		// it gets attributes and uniforms from fragmentCommon3D

		uniform sampler2D tex;

		void main()
		{
			vec4 texel = texture(tex, passTexCoord);

			// apply gamma correction and intensity
			texel.rgb *= intensity;
			outColor.rgb = pow(texel.rgb, vec3(gamma));

			// Apply fog if enabled
			if (fogColor.w > 0.0)
			{
				float depth = gl_FragCoord.z / gl_FragCoord.w;
				float d = fogColor.w * depth;
				float fogFactor = 1.0 - exp(-(d * d));
				outColor.rgb = mix(outColor.rgb, fogColor.rgb, fogFactor);
			}

			outColor.a = texel.a*alpha; // I think alpha shouldn't be modified by gamma and intensity
		}
);

static const char* fragmentSrc3DspriteAlpha = MULTILINE_STRING(

		// it gets attributes and uniforms from fragmentCommon3D

		uniform sampler2D tex;

		void main()
		{
			vec4 texel = texture(tex, passTexCoord);

			if (texel.a <= 0.666)
				discard;

			// apply gamma correction and intensity
			texel.rgb *= intensity;
			outColor.rgb = pow(texel.rgb, vec3(gamma));

			// Apply fog if enabled
			if (fogColor.w > 0.0)
			{
				float depth = gl_FragCoord.z / gl_FragCoord.w;
				float d = fogColor.w * depth;
				float fogFactor = 1.0 - exp(-(d * d));
				outColor.rgb = mix(outColor.rgb, fogColor.rgb, fogFactor);
			}

			outColor.a = texel.a; // in this case alpha from uni3d shouldn't be used
		}
);

static const char* vertexSrc3Dwater = MULTILINE_STRING(

		// it gets attributes and uniforms from vertexCommon3D
		void main()
		{
			vec3 relPos = position - cameraPos;
			vec3 roundedPos = ps1RoundWorld(relPos);

			passTexCoord = texCoord;
			gl_Position = transProjView * transModel * vec4(roundedPos, 1.0);

#if PS1_SCREEN_SPACE_SNAP
			gl_Position.xy = floor(gl_Position.xy * ps1ScreenRes * 0.5)
			                 / (ps1ScreenRes * 0.5);
#endif
		}
);

static const char* vertexSrcAlias = MULTILINE_STRING(

		// it gets attributes and uniforms from vertexCommon3D

		noperspective out vec4 passColor;

		void main()
		{
			vec3 relPos = position - cameraPos;
			vec3 roundedPos = ps1RoundWorld(relPos);

			passColor = vertColor*overbrightbits;
			passTexCoord = texCoord;
			gl_Position = transProjView* transModel * vec4(roundedPos, 1.0);

#if PS1_SCREEN_SPACE_SNAP
			gl_Position.xy = floor(gl_Position.xy * ps1ScreenRes * 0.5)
			                 / (ps1ScreenRes * 0.5);
#endif
		}
);

static const char* fragmentSrcAlias = MULTILINE_STRING(

		// it gets attributes and uniforms from fragmentCommon3D

		uniform sampler2D tex;

		noperspective in vec4 passColor;

		void main()
		{
			vec4 texel = texture(tex, passTexCoord);

			// apply gamma correction and intensity
			texel.rgb *= intensity;
			texel.a *= alpha; // is alpha even used here?
			texel *= min(vec4(1.5), passColor);

			outColor.rgb = pow(texel.rgb, vec3(gamma));

			// Apply fog if enabled
			if (fogColor.w > 0.0)
			{
				float depth = gl_FragCoord.z / gl_FragCoord.w;
				float d = fogColor.w * depth;
				float fogFactor = 1.0 - exp(-(d * d));
				outColor.rgb = mix(outColor.rgb, fogColor.rgb, fogFactor);
			}

			outColor.a = texel.a; // I think alpha shouldn't be modified by gamma and intensity
		}
);

static const char* fragmentSrcAliasColor = MULTILINE_STRING(

		// it gets attributes and uniforms from fragmentCommon3D

		noperspective in vec4 passColor;

		void main()
		{
			vec4 texel = passColor;

			// apply gamma correction and intensity
			texel.a *= alpha; // is alpha even used here?
			outColor.rgb = pow(texel.rgb, vec3(gamma));

			// Apply fog if enabled
			if (fogColor.w > 0.0)
			{
				float depth = gl_FragCoord.z / gl_FragCoord.w;
				float d = fogColor.w * depth;
				float fogFactor = 1.0 - exp(-(d * d));
				outColor.rgb = mix(outColor.rgb, fogColor.rgb, fogFactor);
			}

			outColor.a = texel.a; // I think alpha shouldn't be modified by gamma and intensity
		}
);

static const char* vertexSrcParticles = MULTILINE_STRING(

		// it gets attributes and uniforms from vertexCommon3D

		noperspective out vec4 passColor;

		void main()
		{
			// NOTE: particles are NOT snapped here so they still look smooth;
			// if you want them to wobble too, uncomment the next two lines.
			// vec3 relPos = ps1RoundWorld(position - cameraPos);
			// vec3 relPos = position - cameraPos;
			vec3 relPos = position - cameraPos;

			passColor = vertColor;
			gl_Position = transProjView * transModel * vec4(relPos, 1.0);

			// abusing texCoord for pointSize, pointDist for particles
			float pointDist = texCoord.y*0.1; // with factor 0.1 it looks good.

			gl_PointSize = texCoord.x/pointDist;
		}
);

static const char* fragmentSrcParticles = MULTILINE_STRING(

		// it gets attributes and uniforms from fragmentCommon3D

		noperspective in vec4 passColor;

		void main()
		{
			vec2 offsetFromCenter = 2.0*(gl_PointCoord - vec2(0.5, 0.5));
			float distSquared = dot(offsetFromCenter, offsetFromCenter);
			if (distSquared > 1.0) // this makes sure the particle is round
				discard;

			vec4 texel = passColor;

			// apply gamma correction and intensity
			outColor.rgb = pow(texel.rgb, vec3(gamma));

			// Apply fog if enabled
			if (fogColor.w > 0.0)
			{
				float depth = gl_FragCoord.z / gl_FragCoord.w;
				float d = fogColor.w * depth;
				float fogFactor = 1.0 - exp(-(d * d));
				outColor.rgb = mix(outColor.rgb, fogColor.rgb, fogFactor);
			}

			// fade out towards the edge
			texel.a *= min(1.0, particleFadeFactor*(1.0 - distSquared));

			outColor.a = texel.a;
		}
);

static const char* fragmentSrcParticlesSquare = MULTILINE_STRING(

		// it gets attributes and uniforms from fragmentCommon3D

		noperspective in vec4 passColor;

		void main()
		{
			// so far we didn't use gamma correction for square particles, but this way
			// uniCommon is referenced so hopefully Intels Ivy Bridge HD4000 GPU driver
			// for Windows stops shitting itself
			outColor.rgb = pow(passColor.rgb, vec3(gamma));

			// Apply fog if enabled
			if (fogColor.w > 0.0)
			{
				float depth = gl_FragCoord.z / gl_FragCoord.w;
				float d = fogColor.w * depth;
				float fogFactor = 1.0 - exp(-(d * d));
				outColor.rgb = mix(outColor.rgb, fogColor.rgb, fogFactor);
			}

			outColor.a = passColor.a;
		}
);

static const char* vertexBloomSrcFullScreen = MULTILINE_STRING(

		in vec2 position; // GL4_ATTRIB_POSITION
		in vec2 texCoord; // GL4_ATTRIB_TEXCOORD

		layout (std140) uniform uni2D
		{
			mat4 trans;
		};

		noperspective out vec2 passTexCoord;

		void main()
		{
			gl_Position = trans * vec4(position, 0.0, 1.0);
			passTexCoord = texCoord;
		}
);

static const char* fragmentBloomBright = MULTILINE_STRING(

		noperspective in vec2 passTexCoord;

		uniform sampler2D tex;
		uniform float threshold;

		out vec4 outColor;

		void main()
		{
			vec3 c = texture(tex, passTexCoord).rgb;
			float lum = max(max(c.r, c.g), c.b);

			if (lum > threshold)
			{
				outColor = vec4(c, 1.0);
			}
			else
			{
				outColor = vec4(0.0, 0.0, 0.0, 1.0);
			}
		}
);

static const char* fragmentBloomBlur = MULTILINE_STRING(

		noperspective in vec2 passTexCoord;

		uniform sampler2D tex;
		uniform vec2 dir;

		out vec4 outColor;

		void main()
		{
			vec3 sum = vec3(0.0);

			float w0 = 0.204164;
			float w1 = 0.304005;
			float w2 = 0.093913;
			float w3 = 0.010381;
			float w4 = 0.000489;

			sum += texture(tex, passTexCoord               ).rgb * w0;
			sum += texture(tex, passTexCoord + dir * 1.0   ).rgb * w1;
			sum += texture(tex, passTexCoord - dir * 1.0   ).rgb * w1;
			sum += texture(tex, passTexCoord + dir * 2.0   ).rgb * w2;
			sum += texture(tex, passTexCoord - dir * 2.0   ).rgb * w2;
			sum += texture(tex, passTexCoord + dir * 3.0   ).rgb * w3;
			sum += texture(tex, passTexCoord - dir * 3.0   ).rgb * w3;
			sum += texture(tex, passTexCoord + dir * 4.0   ).rgb * w4;
			sum += texture(tex, passTexCoord - dir * 4.0   ).rgb * w4;

			outColor = vec4(sum, 1.0);
		}
);

#undef MULTILINE_STRING

enum {
	GL4_BINDINGPOINT_UNICOMMON,
	GL4_BINDINGPOINT_UNI2D,
	GL4_BINDINGPOINT_UNI3D,
	GL4_BINDINGPOINT_UNILIGHTS
};

// ============================================================================
// PS1 rounding helpers (CPU side)
// ============================================================================

#define GL4_MAX_PS1_UNIFORMS 32

typedef struct {
	GLuint prog;
	GLint  gridLoc;
	GLint  camLoc;
} gl4PS1Uniforms_t;

static gl4PS1Uniforms_t s_ps1Uniforms[GL4_MAX_PS1_UNIFORMS];
static int   s_numPS1Uniforms = 0;
static float s_ps1GridSize = 1.0f / 1.0f; // default chunkiness
static float s_cameraPos[3] = { 0.0f, 0.0f, 0.0f };

// Called from initShader3D after a program is linked & bound.
static void
registerPS1Uniforms(GLuint prog)
{
	if (s_numPS1Uniforms >= GL4_MAX_PS1_UNIFORMS)
	{
		Com_Printf("WARNING: too many 3D shaders for PS1 uniform cache!\n");
		return;
	}

	GLint gridLoc = glGetUniformLocation(prog, "ps1GridSize");
	GLint camLoc  = glGetUniformLocation(prog, "cameraPos");

	// If neither exists (e.g. fragment-only or special program) just skip.
	if (gridLoc == -1 && camLoc == -1)
	{
		return;
	}

	s_ps1Uniforms[s_numPS1Uniforms].prog    = prog;
	s_ps1Uniforms[s_numPS1Uniforms].gridLoc = gridLoc;
	s_ps1Uniforms[s_numPS1Uniforms].camLoc  = camLoc;
	++s_numPS1Uniforms;

	if (gridLoc != -1)
	{
		glUniform1f(gridLoc, s_ps1GridSize);
	}
	if (camLoc != -1)
	{
		glUniform3f(camLoc, s_cameraPos[0], s_cameraPos[1], s_cameraPos[2]);
	}
}

// Public API: change the PS1 vertex grid size at runtime.
// Smaller values (e.g. 1/64) = chunkier. Larger (e.g. 1/4) = smoother.
void
GL4_SetPS1Grid(float gridSize)
{
	if (gridSize <= 0.0f)
	{
		gridSize = 1.0f / 16.0f;
	}
	s_ps1GridSize = gridSize;

	GLuint prevProg = gl4state.currentShaderProgram;

	int i;
	for (i = 0; i < s_numPS1Uniforms; ++i)
	{
		if (s_ps1Uniforms[i].gridLoc != -1)
		{
			glUseProgram(s_ps1Uniforms[i].prog);
			glUniform1f(s_ps1Uniforms[i].gridLoc, gridSize);
		}
	}

	if (prevProg != 0)
	{
		GL4_UseProgram(prevProg);
	}
}

// Public API: update camera position for precision on large maps.
// Purely optional -- if you never call this, vertices are rounded in
// world space, which is the classic PS1 look anyway.
void
GL4_SetCameraPosition(float x, float y, float z)
{
	if (s_cameraPos[0] == x && s_cameraPos[1] == y && s_cameraPos[2] == z)
	{
		return;
	}

	s_cameraPos[0] = x;
	s_cameraPos[1] = y;
	s_cameraPos[2] = z;

	GLuint prevProg = gl4state.currentShaderProgram;

	int i;
	for (i = 0; i < s_numPS1Uniforms; ++i)
	{
		if (s_ps1Uniforms[i].camLoc != -1)
		{
			glUseProgram(s_ps1Uniforms[i].prog);
			glUniform3f(s_ps1Uniforms[i].camLoc, x, y, z);
		}
	}

	if (prevProg != 0)
	{
		GL4_UseProgram(prevProg);
	}
}

// ============================================================================

static qboolean
initShader2D(gl4ShaderInfo_t* shaderInfo, const char* vertSrc, const char* fragSrc,
	qboolean uniCommonRequired)
{
	GLuint shaders2D[2] = {0};
	GLuint prog = 0;

	if (shaderInfo->shaderProgram != 0)
	{
		Com_Printf("WARNING: calling %s for gl4ShaderInfo_t that already has a shaderProgram!\n",
			__func__);
		glDeleteProgram(shaderInfo->shaderProgram);
	}

	shaderInfo->shaderProgram = 0;
	shaderInfo->uniLmScalesOrTime = -1;
	shaderInfo->uniVblend = -1;

	shaders2D[0] = CompileShader(GL_VERTEX_SHADER, vertSrc, NULL);
	if (shaders2D[0] == 0)
	{
		return false;
	}

	shaders2D[1] = CompileShader(GL_FRAGMENT_SHADER, fragSrc, NULL);
	if (shaders2D[1] == 0)
	{
		glDeleteShader(shaders2D[0]);
		return false;
	}

	prog = CreateShaderProgram(2, shaders2D);

	glDeleteShader(shaders2D[0]);
	glDeleteShader(shaders2D[1]);

	if (prog == 0)
	{
		return false;
	}

	shaderInfo->shaderProgram = prog;
	GL4_UseProgram(prog);

	// Bind the buffer object to the uniform blocks
	GLuint blockIndex = GL_INVALID_INDEX;
	if (uniCommonRequired)
	{
		blockIndex = glGetUniformBlockIndex(prog, "uniCommon");
	}

	if (blockIndex != GL_INVALID_INDEX)
	{
		GLint blockSize;
		glGetActiveUniformBlockiv(prog, blockIndex, GL_UNIFORM_BLOCK_DATA_SIZE, &blockSize);
		if (blockSize != sizeof(gl4state.uniCommonData))
		{
			Com_Printf("WARNING: OpenGL driver disagrees with us about UBO size of 'uniCommon': %i vs %i\n",
					blockSize, (int)sizeof(gl4state.uniCommonData));

			goto err_cleanup;
		}

		glUniformBlockBinding(prog, blockIndex, GL4_BINDINGPOINT_UNICOMMON);
	}
	else if (uniCommonRequired)
	{
		Com_Printf("WARNING: Couldn't find uniform block index 'uniCommon'\n");
		return false;
	}

	blockIndex = glGetUniformBlockIndex(prog, "uni2D");
	if (blockIndex != GL_INVALID_INDEX)
	{
		GLint blockSize;
		glGetActiveUniformBlockiv(prog, blockIndex, GL_UNIFORM_BLOCK_DATA_SIZE, &blockSize);
		if (blockSize != sizeof(gl4state.uni2DData))
		{
			Com_Printf("WARNING: OpenGL driver disagrees with us about UBO size of 'uni2D'\n");
			goto err_cleanup;
		}

		glUniformBlockBinding(prog, blockIndex, GL4_BINDINGPOINT_UNI2D);
	}
	else
	{
		Com_Printf("WARNING: Couldn't find uniform block index 'uni2D'\n");
		goto err_cleanup;
	}

	shaderInfo->uniLmScalesOrTime = glGetUniformLocation(prog, "time");
	if (shaderInfo->uniLmScalesOrTime != -1)
	{
		glUniform1f(shaderInfo->uniLmScalesOrTime, 0.0f);
	}

	shaderInfo->uniVblend = glGetUniformLocation(prog, "v_blend");
	if (shaderInfo->uniVblend != -1)
	{
		glUniform4f(shaderInfo->uniVblend, 0, 0, 0, 0);
	}

	return true;

err_cleanup:

	glDeleteProgram(prog);

	return false;
}

static qboolean
initShader3D(gl4ShaderInfo_t* shaderInfo, const char* vertSrc, const char* fragSrc)
{
	GLuint shaders3D[2] = {0};
	GLuint prog = 0;
	int i=0;

	if (shaderInfo->shaderProgram != 0)
	{
		Com_Printf("WARNING: calling initShader3D for gl4ShaderInfo_t that already has a shaderProgram!\n");
		glDeleteProgram(shaderInfo->shaderProgram);
	}

	shaderInfo->shaderProgram = 0;
	shaderInfo->uniLmScalesOrTime = -1;
	shaderInfo->uniVblend = -1;

	shaders3D[0] = CompileShader(GL_VERTEX_SHADER, vertexCommon3D, vertSrc);
	if (shaders3D[0] == 0)  return false;

	shaders3D[1] = CompileShader(GL_FRAGMENT_SHADER, fragmentCommon3D, fragSrc);
	if (shaders3D[1] == 0)
	{
		glDeleteShader(shaders3D[0]);
		return false;
	}

	prog = CreateShaderProgram(2, shaders3D);

	if (prog == 0)
	{
		goto err_cleanup;
	}

	GL4_UseProgram(prog);

	// Register the PS1/camera uniforms for this program.
	registerPS1Uniforms(prog);

	// Bind the buffer object to the uniform blocks
	GLuint blockIndex = glGetUniformBlockIndex(prog, "uniCommon");
	if (blockIndex != GL_INVALID_INDEX)
	{
		GLint blockSize;
		glGetActiveUniformBlockiv(prog, blockIndex, GL_UNIFORM_BLOCK_DATA_SIZE, &blockSize);
		if (blockSize != sizeof(gl4state.uniCommonData))
		{
			Com_Printf("WARNING: OpenGL driver disagrees with us about UBO size of 'uniCommon'\n");

			goto err_cleanup;
		}

		glUniformBlockBinding(prog, blockIndex, GL4_BINDINGPOINT_UNICOMMON);
	}
	else
	{
		Com_Printf("WARNING: Couldn't find uniform block index 'uniCommon'\n");

		goto err_cleanup;
	}
	blockIndex = glGetUniformBlockIndex(prog, "uni3D");
	if (blockIndex != GL_INVALID_INDEX)
	{
		GLint blockSize;
		glGetActiveUniformBlockiv(prog, blockIndex, GL_UNIFORM_BLOCK_DATA_SIZE, &blockSize);
		if (blockSize != sizeof(gl4state.uni3DData))
		{
			Com_Printf("WARNING: OpenGL driver disagrees with us about UBO size of 'uni3D'\n");
			Com_Printf("         driver says %d, we expect %d\n", blockSize, (int)sizeof(gl4state.uni3DData));

			goto err_cleanup;
		}

		glUniformBlockBinding(prog, blockIndex, GL4_BINDINGPOINT_UNI3D);
	}
	else
	{
		Com_Printf("WARNING: Couldn't find uniform block index 'uni3D'\n");

		goto err_cleanup;
	}
	blockIndex = glGetUniformBlockIndex(prog, "uniLights");
	if (blockIndex != GL_INVALID_INDEX)
	{
		GLint blockSize;
		glGetActiveUniformBlockiv(prog, blockIndex, GL_UNIFORM_BLOCK_DATA_SIZE, &blockSize);
		if (blockSize != sizeof(gl4state.uniLightsData))
		{
			Com_Printf("WARNING: OpenGL driver disagrees with us about UBO size of 'uniLights'\n");
			Com_Printf("         OpenGL says %d, we say %d\n", blockSize, (int)sizeof(gl4state.uniLightsData));

			goto err_cleanup;
		}

		glUniformBlockBinding(prog, blockIndex, GL4_BINDINGPOINT_UNILIGHTS);
	}
	// else: as uniLights is only used in the LM shaders, it's ok if it's missing

	// make sure texture is GL_TEXTURE0
	GLint texLoc = glGetUniformLocation(prog, "tex");
	if (texLoc != -1)
	{
		glUniform1i(texLoc, 0);
	}

	// ..  and the 4 lightmap texture use GL_TEXTURE1..4
	char lmName[10] = "lightmapX";
	for (i=0; i<4; ++i)
	{
		lmName[8] = '0'+i;
		GLint lmLoc = glGetUniformLocation(prog, lmName);
		if (lmLoc != -1)
		{
			glUniform1i(lmLoc, i+1); // lightmap0 belongs to GL_TEXTURE1, lightmap1 to GL_TEXTURE2 etc
		}
	}

	GLint lmScalesLoc = glGetUniformLocation(prog, "lmScales");
	shaderInfo->uniLmScalesOrTime = lmScalesLoc;
	if (lmScalesLoc != -1)
	{
		shaderInfo->lmScales[0] = HMM_Vec4(1.0f, 1.0f, 1.0f, 1.0f);

		for (i=1; i<4; ++i)  shaderInfo->lmScales[i] = HMM_Vec4(0.0f, 0.0f, 0.0f, 0.0f);

		glUniform4fv(lmScalesLoc, 4, shaderInfo->lmScales[0].Elements);
	}

	shaderInfo->shaderProgram = prog;

	// I think the shaders aren't needed anymore once they're linked into the program
	glDeleteShader(shaders3D[0]);
	glDeleteShader(shaders3D[1]);

	return true;

err_cleanup:

	glDeleteShader(shaders3D[0]);
	glDeleteShader(shaders3D[1]);

	if (prog != 0)
	{
		glDeleteProgram(prog);
	}

	return false;
}

static void initUBOs(void)
{
	gl4state.uniCommonData.gamma = 1.0f/vid_gamma->value;
	gl4state.uniCommonData.intensity = gl4_intensity->value;
	gl4state.uniCommonData.intensity2D = gl4_intensity_2D->value;
	gl4state.uniCommonData.color = HMM_Vec4(1, 1, 1, 1);

	glGenBuffers(1, &gl4state.uniCommonUBO);
	glBindBuffer(GL_UNIFORM_BUFFER, gl4state.uniCommonUBO);
	glBindBufferBase(GL_UNIFORM_BUFFER, GL4_BINDINGPOINT_UNICOMMON, gl4state.uniCommonUBO);
	glBufferData(GL_UNIFORM_BUFFER, sizeof(gl4state.uniCommonData), &gl4state.uniCommonData, GL_DYNAMIC_DRAW);

	// the matrix will be set to something more useful later, before being used
	gl4state.uni2DData.transMat4 = HMM_Mat4();

	glGenBuffers(1, &gl4state.uni2DUBO);
	glBindBuffer(GL_UNIFORM_BUFFER, gl4state.uni2DUBO);
	glBindBufferBase(GL_UNIFORM_BUFFER, GL4_BINDINGPOINT_UNI2D, gl4state.uni2DUBO);
	glBufferData(GL_UNIFORM_BUFFER, sizeof(gl4state.uni2DData), &gl4state.uni2DData, GL_DYNAMIC_DRAW);

	// the matrices will be set to something more useful later, before being used
	gl4state.uni3DData.transProjViewMat4 = HMM_Mat4();
	gl4state.uni3DData.transModelMat4 = gl4_identityMat4;
	gl4state.uni3DData.sscroll = 0.0f;
	gl4state.uni3DData.tscroll = 0.0f;
	gl4state.uni3DData.time = 0.0f;
	gl4state.uni3DData.alpha = 1.0f;
	// gl4_overbrightbits 0 means "no scaling" which is equivalent to multiplying with 1
	gl4state.uni3DData.overbrightbits = (gl4_overbrightbits->value <= 0.0f) ? 1.0f : gl4_overbrightbits->value;
	gl4state.uni3DData.particleFadeFactor = gl4_particle_fade_factor->value;
	gl4state.uni3DData.lightScaleForTurb = 1.0f;

	glGenBuffers(1, &gl4state.uni3DUBO);
	glBindBuffer(GL_UNIFORM_BUFFER, gl4state.uni3DUBO);
	glBindBufferBase(GL_UNIFORM_BUFFER, GL4_BINDINGPOINT_UNI3D, gl4state.uni3DUBO);
	glBufferData(GL_UNIFORM_BUFFER, sizeof(gl4state.uni3DData), &gl4state.uni3DData, GL_DYNAMIC_DRAW);

	glGenBuffers(1, &gl4state.uniLightsUBO);
	glBindBuffer(GL_UNIFORM_BUFFER, gl4state.uniLightsUBO);
	glBindBufferBase(GL_UNIFORM_BUFFER, GL4_BINDINGPOINT_UNILIGHTS, gl4state.uniLightsUBO);
	glBufferData(GL_UNIFORM_BUFFER, sizeof(gl4state.uniLightsData), &gl4state.uniLightsData, GL_DYNAMIC_DRAW);

	gl4state.currentUBO = gl4state.uniLightsUBO;
}

static qboolean
createShaders(void)
{
	// reset the PS1 uniform cache before (re)creating programs
	s_numPS1Uniforms = 0;

	if (!initShader2D(&gl4state.si2D, vertexSrc2D, fragmentSrc2D, true))
	{
		Com_Printf("WARNING: Failed to create shader program for textured 2D rendering!\n");
		return false;
	}

	if (!initShader2D(&gl4state.si2Dtinted, vertexSrc2D, fragmentSrc2Dtinted, true))
	{
		Com_Printf("WARNING: Failed to create shader program for tinted 2D rendering!\n");
		return false;
	}

	if (!initShader2D(&gl4state.si2Dcolor, vertexSrc2Dcolor, fragmentSrc2Dcolor, true))
	{
		Com_Printf("WARNING: Failed to create shader program for color-only 2D rendering!\n");
		return false;
	}

	if (!initShader2D(&gl4state.si2DpostProcess, vertexSrc2D, fragmentSrc2Dpostprocess, false))
	{
		Com_Printf("WARNING: Failed to create shader program to render framebuffer object!\n");
		return false;
	}

	if (!initShader2D(&gl4state.si2DpostProcessWater, vertexSrc2D, fragmentSrc2DpostprocessWater, false))
	{
		Com_Printf("WARNING: Failed to create shader program to render framebuffer object under water!\n");
		return false;
	}

	/* bright */
	if (!initShader2D(&gl4state.si2DbloomBright, vertexBloomSrcFullScreen, fragmentBloomBright, false))
	{
		R_Printf(PRINT_ALL, "%s: bright shader failed\n", __func__);
		return false;
	}

	/* blur */
	if (!initShader2D(&gl4state.si2DbloomBlur, vertexBloomSrcFullScreen, fragmentBloomBlur, false))
	{
		R_Printf(PRINT_ALL, "%s: blur shader failed\n", __func__);
		return false;
	}

	const char* lightmappedFrag = (gl4_colorlight->value == 0.0f)
	                               ? fragmentSrc3DlmNoColor : fragmentSrc3Dlm;

	if (!initShader3D(&gl4state.si3Dlm, vertexSrc3Dlm, lightmappedFrag))
	{
		Com_Printf("WARNING: Failed to create shader program for textured 3D rendering with lightmap!\n");
		return false;
	}

	if (!initShader3D(&gl4state.si3Dtrans, vertexSrc3D, fragmentSrc3D))
	{
		Com_Printf("WARNING: Failed to create shader program for rendering translucent 3D things!\n");
		return false;
	}

	if (!initShader3D(&gl4state.si3DcolorOnly, vertexSrc3D, fragmentSrc3Dcolor))
	{
		Com_Printf("WARNING: Failed to create shader program for flat-colored 3D rendering!\n");
		return false;
	}

	if (!initShader3D(&gl4state.si3Dturb, vertexSrc3Dwater, fragmentSrc3Dwater))
	{
		Com_Printf("WARNING: Failed to create shader program for water rendering!\n");
		return false;
	}

	if (!initShader3D(&gl4state.si3DlmFlow, vertexSrc3DlmFlow, lightmappedFrag))
	{
		Com_Printf("WARNING: Failed to create shader program for scrolling textured 3D rendering with lightmap!\n");
		return false;
	}

	if (!initShader3D(&gl4state.si3DtransFlow, vertexSrc3Dflow, fragmentSrc3D))
	{
		Com_Printf("WARNING: Failed to create shader program for scrolling textured translucent 3D rendering!\n");
		return false;
	}

	if (!initShader3D(&gl4state.si3Dsky, vertexSrc3D, fragmentSrc3Dsky))
	{
		Com_Printf("WARNING: Failed to create shader program for sky rendering!\n");
		return false;
	}

	if (!initShader3D(&gl4state.si3Dsprite, vertexSrc3D, fragmentSrc3Dsprite))
	{
		Com_Printf("WARNING: Failed to create shader program for sprite rendering!\n");
		return false;
	}

	if (!initShader3D(&gl4state.si3DspriteAlpha, vertexSrc3D, fragmentSrc3DspriteAlpha))
	{
		Com_Printf("WARNING: Failed to create shader program for alpha-tested sprite rendering!\n");
		return false;
	}

	if (!initShader3D(&gl4state.si3Dalias, vertexSrcAlias, fragmentSrcAlias))
	{
		Com_Printf("WARNING: Failed to create shader program for rendering textured models!\n");
		return false;
	}

	if (!initShader3D(&gl4state.si3DaliasColor, vertexSrcAlias, fragmentSrcAliasColor))
	{
		Com_Printf("WARNING: Failed to create shader program for rendering flat-colored models!\n");
		return false;
	}

	const char* particleFrag = fragmentSrcParticles;
	if (gl4_particle_square->value != 0.0f)
	{
		particleFrag = fragmentSrcParticlesSquare;
	}

	if (!initShader3D(&gl4state.siParticle, vertexSrcParticles, particleFrag))
	{
		Com_Printf("WARNING: Failed to create shader program for rendering particles!\n");
		return false;
	}

	gl4state.currentShaderProgram = 0;

	return true;
}

qboolean
GL4_InitShaders(void)
{
	initUBOs();

	return createShaders();
}

static void deleteShaders(void)
{
	const gl4ShaderInfo_t siZero = {0};
	for (gl4ShaderInfo_t* si = &gl4state.si2D; si <= &gl4state.siParticle; ++si)
	{
		if (si->shaderProgram != 0)
		{
			glDeleteProgram(si->shaderProgram);
		}

		*si = siZero;
	}

	// Invalidate the PS1 uniform cache; programs are gone.
	s_numPS1Uniforms = 0;
}

void
GL4_ShutdownShaders(void)
{
	deleteShaders();

	// let's (ab)use the fact that all 4 UBO handles are consecutive fields
	// of the gl4state struct
	glDeleteBuffers(4, &gl4state.uniCommonUBO);
	gl4state.uniCommonUBO = gl4state.uni2DUBO = gl4state.uni3DUBO = gl4state.uniLightsUBO = 0;
}

qboolean
GL4_RecreateShaders(void)
{
	// delete and recreate the existing shaders (but not the UBOs)
	deleteShaders();
	return createShaders();
}

static inline void
updateUBO(GLuint ubo, GLsizeiptr size, const void *data)
{
	if (gl4state.currentUBO != ubo)
	{
		gl4state.currentUBO = ubo;
		glBindBuffer(GL_UNIFORM_BUFFER, ubo);
	}

	++gl4_numBufferUniforms;

	/*
		atsb: faster in 4.6 and we can use glMapBufferRange to update the entire buffer at once from the beginning.
		we don't need to use glBindBufferRange and we can leave that alone for now.
	*/
	glBufferData(GL_UNIFORM_BUFFER, size, NULL, GL_STREAM_DRAW); // atsb: GL_STREAM_DRAW

	/*
		atsb: we use GL_MAP_WRITE_BIT here to ensure synchronisation between CPU/GPU
		and to prevent any possible data races in cases of sync loss.

		we don't use the persistent mapping feature yet as that would require
		a bit more work with how the buffers are created and mapped.
	*/
	GLvoid* ptr = glMapBufferRange(GL_UNIFORM_BUFFER, 0, size, GL_MAP_WRITE_BIT);
	memcpy(ptr, data, size);
	glUnmapBuffer(GL_UNIFORM_BUFFER);
}

void
GL4_UpdateUBOCommon(void)
{
	updateUBO(gl4state.uniCommonUBO, sizeof(gl4state.uniCommonData), &gl4state.uniCommonData);
}

void
GL4_UpdateUBO2D(void)
{
	updateUBO(gl4state.uni2DUBO, sizeof(gl4state.uni2DData), &gl4state.uni2DData);
}

void
GL4_UpdateUBO3D(void)
{
	updateUBO(gl4state.uni3DUBO, sizeof(gl4state.uni3DData), &gl4state.uni3DData);
}

void
GL4_UpdateUBOLights(void)
{
	updateUBO(gl4state.uniLightsUBO, sizeof(gl4state.uniLightsData), &gl4state.uniLightsData);
}
