/*
 * Copyright (C) 1997-2001 Id Software, Inc.
 * Copyright (C) 2018-2019 Krzysztof Kondrak
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
 * Mesh handling
 *
 * =======================================================================
 */

#include "header/local.h"

#define NUMVERTEXNORMALS 162
#define SHADEDOT_QUANT 16

enum {
	TRIANGLE_STRIP = 0,
	TRIANGLE_FAN = 1
} pipelineIdx;

typedef struct {
	vec3_t vertex;
	float color[4];
	float texCoord[2];
} modelvert;

typedef struct {
	int vertexCount;
	int firstVertex;
} drawinfo_t;

mvtx_t 	*verts_buffer = NULL;
static	drawinfo_t	*drawInfo[2] = {NULL, NULL};
static	modelvert	*vertList[2] = {NULL, NULL};
static	vec3_t	*shadowverts = NULL;
static	int	verts_count = 0;

/* precalculated dot products for quantized angles */
static float r_avertexnormal_dots[SHADEDOT_QUANT][256] = {
#include "../constants/anormtab.h"
};

vec3_t shadevector;
float shadelight[3];
float *shadedots = r_avertexnormal_dots[0];

// correction matrix with "hacked depth" for models with RF_DEPTHHACK flag set
static float r_vulkan_correction_dh[16] = {
	1.f,  0.f, 0.f, 0.f,
	0.f, -1.f, 0.f, 0.f,
	0.f,  0.f, .3f, 0.f,
	0.f,  0.f, .3f, 1.f
};

int
Mesh_VertsRealloc(int count)
{
	void *ptr;

	if (verts_count > count)
	{
		return 0;
	}

	verts_count = ROUNDUP(count * 2, 256);

	ptr = realloc(shadowverts, verts_count * sizeof(vec3_t));
	if (!ptr)
	{
		return -1;
	}
	shadowverts = ptr;

	ptr = realloc(verts_buffer, verts_count * sizeof(mvtx_t));
	if (!ptr)
	{
		return -1;
	}
	verts_buffer = ptr;

	ptr = realloc(vertList[0], verts_count * sizeof(modelvert));
	if (!ptr)
	{
		return -1;
	}
	vertList[0] = ptr;

	ptr = realloc(vertList[1], verts_count * sizeof(modelvert));
	if (!ptr)
	{
		return -1;
	}
	vertList[1] = ptr;

	ptr = realloc(drawInfo[0], verts_count * sizeof(drawinfo_t));
	if (!ptr)
	{
		return -1;
	}
	drawInfo[0] = ptr;

	ptr = realloc(drawInfo[1], verts_count * sizeof(drawinfo_t));
	if (!ptr)
	{
		return -1;
	}
	drawInfo[1] = ptr;

	return 0;
}

/*
===============
Mesh_Init
===============
*/
void
Mesh_Init(void)
{
	shadowverts = NULL;
	verts_buffer = NULL;
	vertList[0] = NULL;
	vertList[1] = NULL;
	drawInfo[0] = NULL;
	drawInfo[1] = NULL;

	verts_count = 0;

	if (Mesh_VertsRealloc(MAX_VERTS))
	{
		Com_Error(ERR_FATAL, "%s: can't allocate memory", __func__);
	}
}

/*
================
Mesh_Free
================
*/
void
Mesh_Free(void)
{
	if (r_validation->value > 1)
	{
		R_Printf(PRINT_ALL, "%s: Deallocated %d mesh verts\n",
			__func__, verts_count);
	}
	verts_count = 0;

	if (shadowverts)
	{
		free(shadowverts);
	}
	shadowverts = NULL;

	if (verts_buffer)
	{
		free(verts_buffer);
	}
	verts_buffer = NULL;

	if (vertList[0])
	{
		free(vertList[0]);
	}
	if (vertList[1])
	{
		free(vertList[1]);
	}
	vertList[0] = NULL;
	vertList[1] = NULL;

	if (drawInfo[0])
	{
		free(drawInfo[0]);
	}
	if (drawInfo[1])
	{
		free(drawInfo[1]);
	}
	drawInfo[0] = NULL;
	drawInfo[1] = NULL;
}

static void
Vk_DrawAliasFrameLerpCommands (entity_t *currententity, int *order, int *order_end,
	float alpha, image_t *skin, float *modelMatrix, int leftHandOffset, int translucentIdx,
	dxtrivertx_t *verts, vec4_t *s_lerped, int verts_count)
{
	int vertCounts[2] = { 0, 0 };
	int pipeCounters[2] = { 0, 0 };
	VkDeviceSize maxTriangleFanIdxCnt = 0;

	if (Mesh_VertsRealloc(verts_count))
	{
		Com_Error(ERR_FATAL, "%s: can't allocate memory", __func__);
	}

	drawInfo[0][0].firstVertex = 0;
	drawInfo[1][0].firstVertex = 0;

	struct {
		float model[16];
		int textured;
	} meshUbo;

	while (1)
	{
		int count;

		/* get the vertex count and primitive type */
		count = *order++;

		if (!count || order >= order_end)
		{
			break; /* done */
		}

		if (count < 0)
		{
			count = -count;
			pipelineIdx = TRIANGLE_FAN;
		}
		else
		{
			pipelineIdx = TRIANGLE_STRIP;
		}

		if (Mesh_VertsRealloc(pipeCounters[pipelineIdx]))
		{
			Com_Error(ERR_FATAL, "%s: can't allocate memory", __func__);
		}

		drawInfo[pipelineIdx][pipeCounters[pipelineIdx]].vertexCount = count;
		maxTriangleFanIdxCnt = Q_max(maxTriangleFanIdxCnt, ((count - 2) * 3));

		if (currententity->flags &
			(RF_SHELL_RED | RF_SHELL_GREEN | RF_SHELL_BLUE))
		{
			meshUbo.textured = 0;
			do
			{
				int vertIdx = vertCounts[pipelineIdx];
				int index_xyz = order[2];

				if (Mesh_VertsRealloc(vertIdx))
				{
					Com_Error(ERR_FATAL, "%s: can't allocate memory", __func__);
				}

				// unused in this case, since texturing is disabled
				vertList[pipelineIdx][vertIdx].texCoord[0] = 0.f;
				vertList[pipelineIdx][vertIdx].texCoord[1] = 0.f;

				vertList[pipelineIdx][vertIdx].color[0] = shadelight[0];
				vertList[pipelineIdx][vertIdx].color[1] = shadelight[1];
				vertList[pipelineIdx][vertIdx].color[2] = shadelight[2];
				vertList[pipelineIdx][vertIdx].color[3] = alpha;

				if (verts_count <= index_xyz)
				{
					R_Printf(PRINT_ALL, "%s: Model has issues with lerped index\n", __func__);
					return;
				}

				vertList[pipelineIdx][vertIdx].vertex[0] = s_lerped[index_xyz][0];
				vertList[pipelineIdx][vertIdx].vertex[1] = s_lerped[index_xyz][1];
				vertList[pipelineIdx][vertIdx].vertex[2] = s_lerped[index_xyz][2];
				vertCounts[pipelineIdx]++;
				order += 3;
			} while (--count);
		}
		else
		{
			/* Do not apply texture for lighmap debug case */
			meshUbo.textured = r_lightmap->value ? 0 : 1;
			do
			{
				int vertIdx = vertCounts[pipelineIdx];
				int index_xyz = order[2];
				float l;

				if (Mesh_VertsRealloc(vertIdx))
				{
					Com_Error(ERR_FATAL, "%s: can't allocate memory", __func__);
				}

				// texture coordinates come from the draw list
				vertList[pipelineIdx][vertIdx].texCoord[0] = ((float *)order)[0];
				vertList[pipelineIdx][vertIdx].texCoord[1] = ((float *)order)[1];

				// normals and vertexes come from the frame list
				l = shadedots[verts[index_xyz].lightnormalindex];

				vertList[pipelineIdx][vertIdx].color[0] = l * shadelight[0];
				vertList[pipelineIdx][vertIdx].color[1] = l * shadelight[1];
				vertList[pipelineIdx][vertIdx].color[2] = l * shadelight[2];
				vertList[pipelineIdx][vertIdx].color[3] = alpha;

				if (verts_count <= index_xyz)
				{
					R_Printf(PRINT_ALL, "%s: Model has issues with lerped index\n", __func__);
					return;
				}

				vertList[pipelineIdx][vertIdx].vertex[0] = s_lerped[index_xyz][0];
				vertList[pipelineIdx][vertIdx].vertex[1] = s_lerped[index_xyz][1];
				vertList[pipelineIdx][vertIdx].vertex[2] = s_lerped[index_xyz][2];
				vertCounts[pipelineIdx]++;
				order += 3;
			} while (--count);
		}

		if (Mesh_VertsRealloc(pipeCounters[pipelineIdx] + 1))
		{
			Com_Error(ERR_FATAL, "%s: can't allocate memory", __func__);
		}

		pipeCounters[pipelineIdx]++;
		drawInfo[pipelineIdx][pipeCounters[pipelineIdx]].firstVertex = vertCounts[pipelineIdx];
	}

	uint32_t uboOffset;
	VkDescriptorSet uboDescriptorSet;
	uint8_t *uboData = QVk_GetUniformBuffer(sizeof(meshUbo), &uboOffset, &uboDescriptorSet);
	memcpy(meshUbo.model, modelMatrix, sizeof(float) * 16);
	memcpy(uboData, &meshUbo, sizeof(meshUbo));

	// player configuration screen model is using the UI renderpass
	int pidx = (r_newrefdef.rdflags & RDF_NOWORLDMODEL) ? RP_UI : RP_WORLD;
	// non-depth write alias models don't occur with RF_WEAPONMODEL set, so no need for additional left-handed pipelines
	qvkpipeline_t pipelines[2][2] = {
		{ vk_drawModelPipelineFan[pidx], vk_drawLefthandModelPipelineFan },
		{ vk_drawNoDepthModelPipelineFan, vk_drawLefthandModelPipelineFan } };
	for (int p = 0; p < 2; p++)
	{
		VkDeviceSize vaoSize = sizeof(modelvert) * vertCounts[p];
		VkBuffer vbo;
		VkDeviceSize vboOffset;
		uint8_t *vertData = QVk_GetVertexBuffer(vaoSize, &vbo, &vboOffset);
		memcpy(vertData, vertList[p], vaoSize);

		QVk_BindPipeline(&pipelines[translucentIdx][leftHandOffset]);
		VkDescriptorSet descriptorSets[] = {
			skin->vk_texture.descriptorSet,
			uboDescriptorSet
		};
		vkCmdBindDescriptorSets(vk_activeCmdbuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, pipelines[translucentIdx][leftHandOffset].layout, 0, 2, descriptorSets, 1, &uboOffset);
		vkCmdBindVertexBuffers(vk_activeCmdbuffer, 0, 1, &vbo, &vboOffset);

		if (p == TRIANGLE_STRIP)
		{
			int i;

			vkCmdBindIndexBuffer(vk_activeCmdbuffer, QVk_GetTriangleStripIbo(maxTriangleFanIdxCnt), 0, VK_INDEX_TYPE_UINT16);

			for (i = 0; i < pipeCounters[p]; i++)
			{
				vkCmdDrawIndexed(vk_activeCmdbuffer, (drawInfo[p][i].vertexCount - 2) * 3, 1, 0, drawInfo[p][i].firstVertex, 0);
			}
		}
		else
		{
			int i;

			vkCmdBindIndexBuffer(vk_activeCmdbuffer, QVk_GetTriangleFanIbo(maxTriangleFanIdxCnt), 0, VK_INDEX_TYPE_UINT16);

			for (i = 0; i < pipeCounters[p]; i++)
			{
				vkCmdDrawIndexed(vk_activeCmdbuffer, (drawInfo[p][i].vertexCount - 2) * 3, 1, 0, drawInfo[p][i].firstVertex, 0);
			}
		}
	}
}

/*
=============
Vk_DrawAliasFrameLerp

interpolates between two frames and origins
FIXME: batch lerp all vertexes
=============
*/
static void
Vk_DrawAliasFrameLerp(entity_t *currententity, dmdx_t *paliashdr, float backlerp, image_t *skin,
	float *modelMatrix, int leftHandOffset, int translucentIdx, vec4_t *s_lerped)
{
	daliasxframe_t *frame, *oldframe;
	dxtrivertx_t *v, *ov, *verts;
	int *order;
	float frontlerp;
	float alpha;
	vec3_t move, delta, vectors[3];
	vec3_t frontv, backv;
	int i;
	int num_mesh_nodes;
	dmdxmesh_t *mesh_nodes;
	qboolean colorOnly = 0 != (currententity->flags &
			(RF_SHELL_RED | RF_SHELL_GREEN | RF_SHELL_BLUE | RF_SHELL_DOUBLE |
			 RF_SHELL_HALF_DAM));

	frame = (daliasxframe_t *)((byte *)paliashdr + paliashdr->ofs_frames
							  + currententity->frame * paliashdr->framesize);
	verts = v = frame->verts;

	oldframe = (daliasxframe_t *)((byte *)paliashdr + paliashdr->ofs_frames
				+ currententity->oldframe * paliashdr->framesize);
	ov = oldframe->verts;

	order = (int *)((byte *)paliashdr + paliashdr->ofs_glcmds);

	if (currententity->flags & RF_TRANSLUCENT)
	{
		alpha = currententity->alpha;
	}
	else
	{
		alpha = 1.0;
	}

	frontlerp = 1.0 - backlerp;

	/* move should be the delta back to the previous frame * backlerp */
	VectorSubtract(currententity->oldorigin, currententity->origin, delta);
	AngleVectors(currententity->angles, vectors[0], vectors[1], vectors[2]);

	move[0] = DotProduct(delta, vectors[0]); /* forward */
	move[1] = -DotProduct(delta, vectors[1]); /* left */
	move[2] = DotProduct(delta, vectors[2]); /* up */

	VectorAdd(move, oldframe->translate, move);

	for (i = 0; i < 3; i++)
	{
		move[i] = backlerp * move[i] + frontlerp * frame->translate[i];

		frontv[i] = frontlerp * frame->scale[i];
		backv[i] = backlerp * oldframe->scale[i];
	}

	R_LerpVerts(colorOnly, paliashdr->num_xyz, v, ov, verts, (float*)s_lerped,
		move, frontv, backv);

	num_mesh_nodes = paliashdr->num_meshes;
	mesh_nodes = (dmdxmesh_t *)((char*)paliashdr + paliashdr->ofs_meshes);

	for (i = 0; i < num_mesh_nodes; i++)
	{
		Vk_DrawAliasFrameLerpCommands(currententity,
			order + mesh_nodes[i].start,
			order + Q_min(paliashdr->num_glcmds,
				mesh_nodes[i].start + mesh_nodes[i].num),
			alpha, skin,
			modelMatrix, leftHandOffset, translucentIdx, verts,
			s_lerped, paliashdr->num_xyz);
	}
}

static void
Vk_DrawAliasShadow(int *order, int *order_end, int posenum,
	float *modelMatrix, entity_t *currententity, vec4_t *s_lerped)
{
	vec3_t	point;
	float	height, lheight;

	lheight = currententity->origin[2] - lightspot[2];

	height = 0;

	height = -lheight + 1.0;

	uint32_t uboOffset;
	VkDescriptorSet uboDescriptorSet;
	uint8_t *uboData = QVk_GetUniformBuffer(sizeof(float) * 16, &uboOffset, &uboDescriptorSet);
	memcpy(uboData, modelMatrix, sizeof(float) * 16);

	while (1)
	{
		int i, count;

		i = 0;
		// get the vertex count and primitive type
		count = *order++;
		if (!count || order >= order_end)
		{
			break; /* done */
		}

		if (count < 0)
		{
			count = -count;
			pipelineIdx = TRIANGLE_FAN;
		}
		else
		{
			pipelineIdx = TRIANGLE_STRIP;
		}

		if (Mesh_VertsRealloc(count))
		{
			Com_Error(ERR_FATAL, "%s: can't allocate memory", __func__);
		}

		do
		{
			if (Mesh_VertsRealloc(order[2]))
			{
				Com_Error(ERR_FATAL, "%s: can't allocate memory", __func__);
			}

			/* normals and vertexes come from the frame list */
			memcpy(point, s_lerped[order[2]], sizeof(point));

			point[0] -= shadevector[0] * (point[2] + lheight);
			point[1] -= shadevector[1] * (point[2] + lheight);
			point[2] = height;

			shadowverts[i][0] = point[0];
			shadowverts[i][1] = point[1];
			shadowverts[i][2] = point[2];

			order += 3;
			i++;
		}
		while (--count);

		if (i > 0)
		{
			VkDeviceSize vaoSize = sizeof(vec3_t) * i;
			VkBuffer vbo;
			VkDeviceSize vboOffset;
			uint8_t *vertData = QVk_GetVertexBuffer(vaoSize, &vbo, &vboOffset);
			memcpy(vertData, shadowverts, vaoSize);

			QVk_BindPipeline(&vk_shadowsPipelineFan);
			vkCmdBindDescriptorSets(vk_activeCmdbuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, vk_shadowsPipelineFan.layout, 0, 1, &uboDescriptorSet, 1, &uboOffset);
			vkCmdBindVertexBuffers(vk_activeCmdbuffer, 0, 1, &vbo, &vboOffset);

			if (pipelineIdx == TRIANGLE_STRIP)
			{
				vkCmdBindIndexBuffer(vk_activeCmdbuffer, QVk_GetTriangleStripIbo((i - 2) * 3), 0, VK_INDEX_TYPE_UINT16);
				vkCmdDrawIndexed(vk_activeCmdbuffer, (i - 2) * 3, 1, 0, 0, 0);
			}
			else
			{
				vkCmdBindIndexBuffer(vk_activeCmdbuffer, QVk_GetTriangleFanIbo((i - 2) * 3), 0, VK_INDEX_TYPE_UINT16);
				vkCmdDrawIndexed(vk_activeCmdbuffer, (i - 2) * 3, 1, 0, 0, 0);
			}
		}
	}
}

static qboolean
R_CullAliasModel(const model_t *currentmodel, vec3_t bbox[8], entity_t *e)
{
	dmdx_t *paliashdr;

	paliashdr = (dmdx_t *)currentmodel->extradata;
	if (!paliashdr)
	{
		R_Printf(PRINT_ALL, "%s %s: Model is not fully loaded\n",
				__func__, currentmodel->name);
		return true;
	}

	if ((e->frame >= paliashdr->num_frames) || (e->frame < 0))
	{
		R_Printf(PRINT_DEVELOPER, "%s %s: no such frame %d\n",
				__func__, currentmodel->name, e->frame);
		e->frame = 0;
	}

	if ((e->oldframe >= paliashdr->num_frames) || (e->oldframe < 0))
	{
		R_Printf(PRINT_DEVELOPER, "%s %s: no such oldframe %d\n",
				__func__, currentmodel->name, e->oldframe);
		e->oldframe = 0;
	}

	return R_CullAliasMeshModel(paliashdr, frustum, e->frame, e->oldframe,
		e->angles, e->origin, bbox);
}

void
R_DrawAliasModel(entity_t *currententity, const model_t *currentmodel)
{
	int leftHandOffset = 0, i;
	float prev_viewproj[16], an;
	dmdx_t *paliashdr;
	vec4_t *s_lerped;

	if (!(currententity->flags & RF_WEAPONMODEL))
	{
		vec3_t bbox[8];

		if (R_CullAliasModel(currentmodel, bbox, currententity))
		{
			return;
		}
	}

	if (currententity->flags & RF_WEAPONMODEL)
	{
		if (r_lefthand->value == 2)
		{
			return;
		}
	}

	paliashdr = (dmdx_t *)currentmodel->extradata;

	/* get lighting information */
	if (currententity->flags &
		(RF_SHELL_HALF_DAM | RF_SHELL_GREEN | RF_SHELL_RED |
		 RF_SHELL_BLUE | RF_SHELL_DOUBLE))
	{
		VectorClear(shadelight);

		if (currententity->flags & RF_SHELL_HALF_DAM)
		{
			shadelight[0] = 0.56;
			shadelight[1] = 0.59;
			shadelight[2] = 0.45;
		}

		if (currententity->flags & RF_SHELL_DOUBLE)
		{
			shadelight[0] = 0.9;
			shadelight[1] = 0.7;
		}

		if (currententity->flags & RF_SHELL_RED)
		{
			shadelight[0] = 1.0;
		}

		if (currententity->flags & RF_SHELL_GREEN)
		{
			shadelight[1] = 1.0;
		}

		if (currententity->flags & RF_SHELL_BLUE)
		{
			shadelight[2] = 1.0;
		}
	}
	else if (currententity->flags & RF_FULLBRIGHT)
	{
		for (i = 0; i < 3; i++)
		{
			shadelight[i] = 1.0;
		}
	}
	else
	{
		if (!r_worldmodel || !r_worldmodel->lightdata)
		{
			shadelight[0] = shadelight[1] = shadelight[2] = 1.0F;
		}
		else
		{
			R_LightPoint(r_worldmodel->grid, currententity, &r_newrefdef, r_worldmodel->surfaces,
				r_worldmodel->nodes, currententity->origin, shadelight,
				r_modulate->value, lightspot);
		}

		/* player lighting hack for communication back to server */
		if (currententity->flags & RF_WEAPONMODEL)
		{
			/* pick the greatest component, which should be
			   the same as the mono value returned by software */
			if (shadelight[0] > shadelight[1])
			{
				if (shadelight[0] > shadelight[2])
				{
					r_lightlevel->value = 150 * shadelight[0];
				}
				else
				{
					r_lightlevel->value = 150 * shadelight[2];
				}
			}
			else
			{
				if (shadelight[1] > shadelight[2])
				{
					r_lightlevel->value = 150 * shadelight[1];
				}
				else
				{
					r_lightlevel->value = 150 * shadelight[2];
				}
			}
		}
	}

	if (currententity->flags & RF_MINLIGHT)
	{
		for (i = 0; i < 3; i++)
		{
			if (shadelight[i] > 0.1)
			{
				break;
			}
		}

		if (i == 3)
		{
			shadelight[0] = 0.1;
			shadelight[1] = 0.1;
			shadelight[2] = 0.1;
		}
	}

	if (currententity->flags & RF_GLOW)
	{
		/* bonus items will pulse with time */
		float scale;

		scale = 0.1 * sin(r_newrefdef.time * 7);

		for (i = 0; i < 3; i++)
		{
			float	min;

			min = shadelight[i] * 0.8;
			shadelight[i] += scale;

			if (shadelight[i] < min)
			{
				shadelight[i] = min;
			}
		}
	}

	/* ir goggles color override */
	if (r_newrefdef.rdflags & RDF_IRGOGGLES && currententity->flags &
		RF_IR_VISIBLE)
	{
		shadelight[0] = 1.0;
		shadelight[1] = 0.0;
		shadelight[2] = 0.0;
	}

	shadedots = r_avertexnormal_dots[((int)(currententity->angles[1] *
				(SHADEDOT_QUANT / 360.0))) & (SHADEDOT_QUANT - 1)];

	an = currententity->angles[1] / 180 * M_PI;
	shadevector[0] = cos(-an);
	shadevector[1] = sin(-an);
	shadevector[2] = 1;
	VectorNormalize(shadevector);

	/* locate the proper data */
	c_alias_polys += paliashdr->num_tris;

	/* draw all the triangles */
	if (currententity->flags & RF_DEPTHHACK || r_newrefdef.rdflags & RDF_NOWORLDMODEL) { // hack the depth range to prevent view model from poking into walls
		float r_proj_aspect = (float)r_newrefdef.width / r_newrefdef.height;
		float r_proj_fovy = r_newrefdef.fov_y;
		float dist = (r_farsee->value == 0) ? 4096.0f : 8192.0f;
		// use different range for player setup screen so it doesn't collide with the viewmodel
		r_vulkan_correction_dh[10] = 0.3f - (r_newrefdef.rdflags & RDF_NOWORLDMODEL) * 0.1f;
		r_vulkan_correction_dh[14] = 0.3f - (r_newrefdef.rdflags & RDF_NOWORLDMODEL) * 0.1f;

		memcpy(prev_viewproj, r_viewproj_matrix, sizeof(r_viewproj_matrix));
		if (currententity->flags & RF_WEAPONMODEL && r_gunfov->value < 0)
		{
			Mat_Perspective(r_projection_matrix, r_vulkan_correction_dh, r_proj_fovy, r_proj_aspect, 4, dist);
		}
		else
		{
			Mat_Perspective(r_projection_matrix, r_vulkan_correction_dh, r_gunfov->value, r_proj_aspect, 4, dist);
		}
		Mat_Mul(r_view_matrix, r_projection_matrix, r_viewproj_matrix);
		vkCmdPushConstants(vk_activeCmdbuffer, vk_drawTexQuadPipeline[vk_state.current_renderpass].layout, VK_SHADER_STAGE_VERTEX_BIT, 0, sizeof(r_viewproj_matrix), r_viewproj_matrix);
	}

	if ( ( currententity->flags & RF_WEAPONMODEL ) && ( r_lefthand->value == 1.0F ) )
	{
		Mat_Scale(r_viewproj_matrix, -1.f, 1.f, 1.f);
		vkCmdPushConstants(vk_activeCmdbuffer, vk_drawTexQuadPipeline[vk_state.current_renderpass].layout, VK_SHADER_STAGE_VERTEX_BIT, 0, sizeof(r_viewproj_matrix), r_viewproj_matrix);
		leftHandOffset = 1;
	}

	currententity->angles[PITCH] = -currententity->angles[PITCH];	// sigh.
	/* buffer for scalled vert from frame */
	s_lerped = R_VertBufferRealloc(paliashdr->num_xyz);

	{
		float model[16];
		image_t *skin = NULL;
		Mat_Identity(model);
		R_RotateForEntity (currententity, model);

		currententity->angles[PITCH] = -currententity->angles[PITCH];	// sigh.

		// select skin
		if (currententity->skin)
		{
			skin = currententity->skin;	// custom player skin
		}
		else
		{
			if (currententity->skinnum < currentmodel->numskins)
			{
				skin = currentmodel->skins[currententity->skinnum];
			}

			if (!skin)
			{
				skin = currentmodel->skins[0];
			}
		}

		if (!skin)
		{
			skin = r_notexture;	// fallback...
		}

		// draw it
		if ( (currententity->frame >= paliashdr->num_frames)
			|| (currententity->frame < 0) )
		{
			R_Printf(PRINT_ALL, "%s %s: no such frame %d\n",
				__func__, currentmodel->name, currententity->frame);
			currententity->frame = 0;
			currententity->oldframe = 0;
		}

		if ( (currententity->oldframe >= paliashdr->num_frames)
			|| (currententity->oldframe < 0))
		{
			R_Printf(PRINT_ALL, "%s %s: no such oldframe %d\n",
				__func__, currentmodel->name, currententity->oldframe);
			currententity->frame = 0;
			currententity->oldframe = 0;
		}

		if ( !r_lerpmodels->value )
		{
			currententity->backlerp = 0;
		}

		Vk_DrawAliasFrameLerp(currententity, paliashdr, currententity->backlerp,
			skin, model, leftHandOffset, (currententity->flags & RF_TRANSLUCENT) ? 1 : 0,
			s_lerped);
	}

	if ( ( currententity->flags & RF_WEAPONMODEL ) && ( r_lefthand->value == 1.0F ) )
	{
		Mat_Scale(r_viewproj_matrix, -1.f, 1.f, 1.f);
		vkCmdPushConstants(vk_activeCmdbuffer, vk_drawTexQuadPipeline[vk_state.current_renderpass].layout, VK_SHADER_STAGE_VERTEX_BIT, 0, sizeof(r_viewproj_matrix), r_viewproj_matrix);
	}

	if (currententity->flags & RF_DEPTHHACK || r_newrefdef.rdflags & RDF_NOWORLDMODEL)
	{
		memcpy(r_viewproj_matrix, prev_viewproj, sizeof(prev_viewproj));
		vkCmdPushConstants(vk_activeCmdbuffer, vk_drawTexQuadPipeline[vk_state.current_renderpass].layout, VK_SHADER_STAGE_VERTEX_BIT, 0, sizeof(r_viewproj_matrix), r_viewproj_matrix);
	}

	if (vk_shadows->value && !(currententity->flags & (RF_TRANSLUCENT | RF_WEAPONMODEL)))
	{
		int num_mesh_nodes, i;
		dmdxmesh_t *mesh_nodes;
		float model[16];
		int *order;

		Mat_Identity(model);
		R_RotateForEntity(currententity, model);

		order = (int *)((byte *)paliashdr + paliashdr->ofs_glcmds);

		num_mesh_nodes = paliashdr->num_meshes;
		mesh_nodes = (dmdxmesh_t *)((char*)paliashdr + paliashdr->ofs_meshes);

		for (i = 0; i < num_mesh_nodes; i++)
		{
			Vk_DrawAliasShadow (
				order + mesh_nodes[i].start,
				order + Q_min(paliashdr->num_glcmds,
					mesh_nodes[i].start + mesh_nodes[i].num),
				currententity->frame, model, currententity,
				s_lerped);
		}
	}
}
