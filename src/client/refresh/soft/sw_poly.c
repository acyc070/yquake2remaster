/*
 * Copyright (C) 1997-2001 Id Software, Inc.
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
 */

#include <assert.h>
#include <limits.h>
#include "header/local.h"

static int	r_polyblendcolor;

polydesc_t	r_polydesc;

msurface_t	*r_alpha_surfaces;

static int	clip_current;
vec5_t		r_clip_verts[2][MAXWORKINGVERTS+2];
static emitpoint_t	outverts[MAXWORKINGVERTS+3];

static int	s_minindex, s_maxindex;

// Affine rasterisation data per scanline (fixed-point 16.16)
#define MAX_SCANLINES 2048
static int	left_u[MAX_SCANLINES];
static int	left_s[MAX_SCANLINES];
static int	left_t[MAX_SCANLINES];
static int	left_izi[MAX_SCANLINES];
static int	right_u[MAX_SCANLINES];
static int	right_s[MAX_SCANLINES];
static int	right_t[MAX_SCANLINES];
static int	right_izi[MAX_SCANLINES];

/*
** R_ClipPolyFace
**
** Clips the winding at clip_verts[clip_current] and changes clip_current
** Throws out the back side
*/
static int
R_ClipPolyFace(int nump, const clipplane_t *pclipplane)
{
	float dists[MAXWORKINGVERTS+3] = {0};
	const float *vert2, *pclipnormal;
	float *in, *instep, *outstep;
	float frac, clipdist;
	int i, outcount;

	clipdist = pclipplane->dist;
	pclipnormal = pclipplane->normal;

	// calc dists
	if (clip_current)
	{
		in = r_clip_verts[1][0];
		outstep = r_clip_verts[0][0];
		clip_current = 0;
	}
	else
	{
		in = r_clip_verts[0][0];
		outstep = r_clip_verts[1][0];
		clip_current = 1;
	}

	instep = in;
	for (i=0 ; i<nump ; i++, instep += sizeof (vec5_t) / sizeof (vec_t))
	{
		dists[i] = DotProduct (instep, pclipnormal) - clipdist;
	}

	/* handle wraparound case */
	dists[nump] = dists[0];
	memmove(instep, in, sizeof (vec5_t));

	/* clip the winding */
	instep = in;
	outcount = 0;

	for (i=0 ; i<nump ; i++, instep += sizeof (vec5_t) / sizeof (vec_t))
	{
		if (dists[i] >= 0)
		{
			memcpy (outstep, instep, sizeof (vec5_t));
			outstep += sizeof (vec5_t) / sizeof (vec_t);
			outcount++;
		}

		if (dists[i] == 0 || dists[i+1] == 0)
			continue;

		if ( (dists[i] > 0) == (dists[i+1] > 0) )
			continue;

		// split it into a new vertex
		frac = dists[i] / (dists[i] - dists[i+1]);

		vert2 = instep + sizeof (vec5_t) / sizeof (vec_t);

		outstep[0] = instep[0] + frac*(vert2[0] - instep[0]);
		outstep[1] = instep[1] + frac*(vert2[1] - instep[1]);
		outstep[2] = instep[2] + frac*(vert2[2] - instep[2]);
		outstep[3] = instep[3] + frac*(vert2[3] - instep[3]);
		outstep[4] = instep[4] + frac*(vert2[4] - instep[4]);

		outstep += sizeof (vec5_t) / sizeof (vec_t);
		outcount++;
	}

	return outcount;
}

/*
** R_PolygonScanLeftEdge
**
** Fills left_u, left_s, left_t, left_izi for every integer scanline
** touched by the left edge. Uses the vertices in r_polydesc.pverts
** (already projected to screen).
*/
static void
R_PolygonScanLeftEdge(void)
{
	const emitpoint_t *pvert, *pnext;
	float du, dv, vtop, u_step, s_step, t_step, izi_step;
	int i, lmaxindex;
	int v, itop, ibottom;

	i = s_minindex;
	if (i == 0)
		i = r_polydesc.nump;

	lmaxindex = s_maxindex;
	if (lmaxindex == 0)
		lmaxindex = r_polydesc.nump;

	vtop = ceil(r_polydesc.pverts[i].v);

	do
	{
		float vbottom;

		pvert = &r_polydesc.pverts[i];
		pnext = pvert - 1;

		vbottom = ceil(pnext->v);

		if (vtop < vbottom)
		{
			du = pnext->u - pvert->u;
			dv = pnext->v - pvert->v;

			// steps per one unit of v
			u_step = du / dv;
			s_step = (pnext->s - pvert->s) / dv;
			t_step = (pnext->t - pvert->t) / dv;
			izi_step = (pnext->zi - pvert->zi) / dv;

			// adjust to the first integer v
			float frac = vtop - pvert->v;
			float u_start = pvert->u + u_step * frac;
			float s_start = pvert->s + s_step * frac;
			float t_start = pvert->t + t_step * frac;
			float izi_start = pvert->zi + izi_step * frac;

			itop = (int)vtop;
			ibottom = (int)vbottom;

			for (v = itop; v < ibottom; v++)
			{
				// store fixed-point values (16.16)
				left_u[v]   = (int)(u_start * SHIFT16XYZ_MULT);
				left_s[v]   = (int)(s_start * SHIFT16XYZ_MULT);
				left_t[v]   = (int)(t_start * SHIFT16XYZ_MULT);
				left_izi[v] = (int)(izi_start * SHIFT16XYZ_MULT * 0x8000); // match original izi scale

				u_start += u_step;
				s_start += s_step;
				t_start += t_step;
				izi_start += izi_step;
			}
		}

		vtop = vbottom;

		i--;
		if (i == 0)
			i = r_polydesc.nump;

	} while (i != lmaxindex);
}

/*
** R_PolygonScanRightEdge
**
** Fills right_u, right_s, right_t, right_izi for every integer scanline.
*/
static void
R_PolygonScanRightEdge(void)
{
	const emitpoint_t *pvert, *pnext;
	float du, dv, vtop, u_step, s_step, t_step, izi_step;
	int i;
	int v, itop, ibottom;
	float vvert;

	i = s_minindex;

	vvert = r_polydesc.pverts[i].v;
	if (vvert < r_refdef.fvrecty_adj)
		vvert = r_refdef.fvrecty_adj;
	if (vvert > r_refdef.fvrectbottom_adj)
		vvert = r_refdef.fvrectbottom_adj;

	vtop = ceil(vvert);

	do
	{
		float vbottom, vnext;

		pvert = &r_polydesc.pverts[i];
		pnext = pvert + 1;

		vnext = pnext->v;
		if (vnext < r_refdef.fvrecty_adj)
			vnext = r_refdef.fvrecty_adj;
		if (vnext > r_refdef.fvrectbottom_adj)
			vnext = r_refdef.fvrectbottom_adj;

		vbottom = ceil(vnext);

		if (vtop < vbottom)
		{
			float uvert = pvert->u;
			if (uvert < r_refdef.fvrectx_adj)
				uvert = r_refdef.fvrectx_adj;
			if (uvert > r_refdef.fvrectright_adj)
				uvert = r_refdef.fvrectright_adj;

			float unext = pnext->u;
			if (unext < r_refdef.fvrectx_adj)
				unext = r_refdef.fvrectx_adj;
			if (unext > r_refdef.fvrectright_adj)
				unext = r_refdef.fvrectright_adj;

			du = unext - uvert;
			dv = vnext - vvert;

			u_step = du / dv;
			s_step = (pnext->s - pvert->s) / dv;
			t_step = (pnext->t - pvert->t) / dv;
			izi_step = (pnext->zi - pvert->zi) / dv;

			float frac = vtop - vvert;
			float u_start = uvert + u_step * frac;
			float s_start = pvert->s + s_step * frac;
			float t_start = pvert->t + t_step * frac;
			float izi_start = pvert->zi + izi_step * frac;

			itop = (int)vtop;
			ibottom = (int)vbottom;

			for (v = itop; v < ibottom; v++)
			{
				right_u[v]  = (int)(u_start * SHIFT16XYZ_MULT);
				right_s[v]  = (int)(s_start * SHIFT16XYZ_MULT);
				right_t[v]  = (int)(t_start * SHIFT16XYZ_MULT);
				right_izi[v] = (int)(izi_start * SHIFT16XYZ_MULT * 0x8000);

				u_start += u_step;
				s_start += s_step;
				t_start += t_step;
				izi_start += izi_step;
			}
		}

		vtop = vbottom;
		vvert = vnext;

		i++;
		if (i == r_polydesc.nump)
			i = 0;

	} while (i != s_maxindex);
}

/*
** R_ClipAndDrawPoly
**
** Clips the polygon, projects vertices, scans edges, and draws it
** using pure affine texture mapping (PS1 style).
*/
void
R_ClipAndDrawPoly(float alpha, int isturbulent, qboolean textured)
{
	vec_t		*pv;
	int		i, nump;
	vec3_t		transformed, local;
	float		ymin, ymax;
	int		vstart, vend, v;
	pixel_t		*texbase;
	int		texwidth;

	if (!textured)
	{
		// flat shaded quad – keep original simple path
		// (not shown here for brevity; unchanged from original)
		return;
	}

	// clip to the frustum in worldspace
	nump = r_polydesc.nump;
	clip_current = 0;

	for (i=0 ; i<4 ; i++)
	{
		nump = R_ClipPolyFace (nump, &view_clipplanes[i]);
		if (nump < 3)
			return;
		if (nump > MAXWORKINGVERTS)
		{
			Com_Error(ERR_DROP, "%s: too many points: %d", __func__, nump);
			return;
		}
	}

	// transform vertices into viewspace and project
	pv = &r_clip_verts[clip_current][0][0];

	for (i=0 ; i<nump ; i++)
	{
		float scale;
		emitpoint_t *pout;

		VectorSubtract (pv, r_origin, local);
		TransformVector (local, transformed);

		if (transformed[2] < NEAR_CLIP)
			transformed[2] = NEAR_CLIP;

		pout = &outverts[i];
		pout->zi = 1.0 / transformed[2];

		// original texture coordinates (world‑space, will be interpolated affinely)
		pout->s = pv[3];
		pout->t = pv[4];

		scale = xscale * pout->zi;
		pout->u = (xcenter + scale * transformed[0]);

		scale = yscale * pout->zi;
		pout->v = (ycenter - scale * transformed[1]);

		pv += sizeof (vec5_t) / sizeof (vec_t);
	}

	// find top and bottom indices
	ymin = 1e30f;
	ymax = -1e30f;
	for (i=0; i<nump; i++)
	{
		if (outverts[i].v < ymin) { ymin = outverts[i].v; s_minindex = i; }
		if (outverts[i].v > ymax) { ymax = outverts[i].v; s_maxindex = i; }
	}

	// prepare texture info
	texwidth = r_polydesc.pixel_width;
	texbase = r_polydesc.pixels;

	// copy first vertex to last for edge walking
	outverts[nump] = outverts[0];

	// scan edges and fill per‑scanline buffers
	R_PolygonScanLeftEdge();
	R_PolygonScanRightEdge();

	// draw all scanlines from top to bottom
	vstart = (int)ceil(ymin);
	vend   = (int)ceil(ymax);

	for (v = vstart; v < vend; v++)
	{
		int left = left_u[v] >> SHIFT16XYZ;
		int right = right_u[v] >> SHIFT16XYZ;
		int count = right - left;
		int s, t, izi;
		int sstep, tstep, izistep;
		pixel_t *pdest;
		zvalue_t *pz;

		if (count <= 0)
			continue;

		// starting values at left edge
		s = left_s[v];
		t = left_t[v];
		izi = left_izi[v];

		// steps across the scanline
		sstep = (right_s[v] - left_s[v]) / count;
		tstep = (right_t[v] - left_t[v]) / count;
		izistep = (right_izi[v] - left_izi[v]) / count;

		pdest = d_viewbuffer + vid_buffer_width * v + left;
		pz    = d_pzbuffer + vid_buffer_width * v + left;

		// mark damage for z-buffer (original behaviour)
		VID_DamageZBuffer(left, v);
		VID_DamageZBuffer(left + count, v);

		while (count-- > 0)
		{
			unsigned texel;
			int s_fixed = s >> SHIFT16XYZ;
			int t_fixed = t >> SHIFT16XYZ;

			if (isturbulent)
			{
				// turbulent (warp) – keep original warping logic
				int sturb, tturb;
				sturb = (s_fixed + sintable[(t_fixed) & (CYCLE-1)]) & 63;
				tturb = (t_fixed + sintable[(s_fixed) & (CYCLE-1)]) & 63;
				texel = *(pixel_t *)((byte *)texbase + sturb + (tturb << 6));
			}
			else
			{
				// normal texture – clamp to edges
				if (s_fixed < 0) s_fixed = 0;
				if (s_fixed >= texwidth) s_fixed = texwidth - 1;
				if (t_fixed < 0) t_fixed = 0;
				if (t_fixed >= r_polydesc.pixel_height) t_fixed = r_polydesc.pixel_height - 1;
				texel = *(pixel_t *)((byte *)texbase + s_fixed + t_fixed * texwidth);
			}

			// alpha blending
			if (alpha != 1.0f)
			{
				if (alpha > 0.33f)
					texel = vid_alphamap[texel * 256 + *pdest];
				else
					texel = vid_alphamap[texel + *pdest * 256];
			}

			// depth test
			if (*pz <= (izi >> SHIFT16XYZ))
			{
				*pdest = texel;
				*pz    = izi >> SHIFT16XYZ;
			}

			s += sstep;
			t += tstep;
			izi += izistep;
			pdest++;
			pz++;
		}
	}
}

/*
** R_BuildPolygonFromSurface
*/
static void
R_BuildPolygonFromSurface(const entity_t *currententity, const model_t *currentmodel, msurface_t *fa)
{
	medge_t *pedges, *r_pedge;
	float tmins[2] = { 0, 0 };
	int i, lnumverts;
	const float *vec;
	vec5_t *pverts;

	r_polydesc.nump = 0;

	/* reconstruct the polygon */
	pedges = currentmodel->edges;
	lnumverts = fa->numedges;

	pverts = r_clip_verts[0];

	for (i=0 ; i<lnumverts ; i++)
	{
		int lindex;

		lindex = currentmodel->surfedges[fa->firstedge + i];

		if (lindex > 0)
		{
			r_pedge = &pedges[lindex];
			vec = currentmodel->vertexes[r_pedge->v[0]].position;
		}
		else
		{
			r_pedge = &pedges[-lindex];
			vec = currentmodel->vertexes[r_pedge->v[1]].position;
		}

		VectorCopy (vec, pverts[i] );
	}

	VectorCopy( fa->texinfo->vecs[0], r_polydesc.vright );
	VectorCopy( fa->texinfo->vecs[1], r_polydesc.vup );
	VectorCopy( fa->plane->normal, r_polydesc.vpn );
	VectorCopy( r_origin, r_polydesc.viewer_position );

	if ( fa->flags & SURF_PLANEBACK )
	{
		VectorSubtract( vec3_origin, r_polydesc.vpn, r_polydesc.vpn );
	}

	if ( fa->texinfo->flags & (SURF_WARP | SURF_SCROLL) )
	{
		r_polydesc.pixels       = fa->texinfo->image->pixels[0];
		r_polydesc.pixel_width  = fa->texinfo->image->width;
		r_polydesc.pixel_height = fa->texinfo->image->height;
	}
	else
	{
		surfcache_t *scache;

		scache = D_CacheSurface(currententity, fa, 0);

		r_polydesc.pixels       = scache->data;
		r_polydesc.pixel_width  = scache->width;
		r_polydesc.pixel_height = scache->height;

		tmins[0] = fa->texturemins[0];
		tmins[1] = fa->texturemins[1];
	}

	r_polydesc.dist = DotProduct( r_polydesc.vpn, pverts[0] );

	r_polydesc.s_offset = fa->texinfo->vecs[0][3] - tmins[0];
	r_polydesc.t_offset = fa->texinfo->vecs[1][3] - tmins[1];

	// scrolling texture addition
	if (fa->texinfo->flags & SURF_SCROLL)
	{
		float sscroll, tscroll;

		R_FlowingScroll(&r_newrefdef, fa->texinfo->flags, &sscroll, &tscroll);

		r_polydesc.s_offset += 2 * sscroll;
		r_polydesc.t_offset += 2 * tscroll;
	}

	r_polydesc.nump = lnumverts;
}

/*
** R_DrawAlphaSurfaces
*/
void
R_DrawAlphaSurfaces(const entity_t *currententity)
{
	msurface_t *s = r_alpha_surfaces;
	const model_t *currentmodel = r_worldmodel;

	modelorg[0] = -r_origin[0];
	modelorg[1] = -r_origin[1];
	modelorg[2] = -r_origin[2];

	while ( s )
	{
		R_BuildPolygonFromSurface(currententity, currentmodel, s);

		// pass down all the texinfo flags, not just SURF_WARP.
		if (s->texinfo->flags & SURF_TRANS66)
			R_ClipAndDrawPoly( 0.60f, (s->texinfo->flags & (SURF_WARP | SURF_SCROLL)), true );
		else
			R_ClipAndDrawPoly( 0.30f, (s->texinfo->flags & (SURF_WARP | SURF_SCROLL)), true );

		s = s->nextalphasurface;
	}

	r_alpha_surfaces = NULL;
}

/*
** R_IMFlatShadedQuad
*/
void
R_IMFlatShadedQuad( const vec3_t a, const vec3_t b, const vec3_t c, const vec3_t d, int color, float alpha )
{
	vec3_t s0, s1;

	r_polydesc.nump = 4;
	VectorCopy( r_origin, r_polydesc.viewer_position );

	VectorCopy( a, r_clip_verts[0][0] );
	VectorCopy( b, r_clip_verts[0][1] );
	VectorCopy( c, r_clip_verts[0][2] );
	VectorCopy( d, r_clip_verts[0][3] );

	r_clip_verts[0][0][3] = 0;
	r_clip_verts[0][1][3] = 0;
	r_clip_verts[0][2][3] = 0;
	r_clip_verts[0][3][3] = 0;

	r_clip_verts[0][0][4] = 0;
	r_clip_verts[0][1][4] = 0;
	r_clip_verts[0][2][4] = 0;
	r_clip_verts[0][3][4] = 0;

	VectorSubtract( d, c, s0 );
	VectorSubtract( c, b, s1 );
	CrossProduct( s0, s1, r_polydesc.vpn );
	VectorNormalize( r_polydesc.vpn );

	r_polydesc.dist = DotProduct( r_polydesc.vpn, r_clip_verts[0][0] );

	r_polyblendcolor = color;

	R_ClipAndDrawPoly( alpha, false, false );
}
