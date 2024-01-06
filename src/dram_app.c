#ifdef DEBUG
 	#include <stdio.h>
#endif
#ifdef MAIN
 	#include <stdio.h>
 	#include <string.h>
#endif

#include <stdlib.h>						
#include "../include/sea.h"				
#include "../util/util.h"						


struct mpos
DECLARE_FUNC(diag_affine_banded_fill, SUFFIX)(
	struct sea_result *aln,
	struct sea_params param,
	char *mat);

struct mpos
DECLARE_FUNC(diag_affine_banded_search, SUFFIX)(
	struct sea_result *aln,
	struct sea_params param,
	char *mat,
	struct mpos o);

sea_int_t
DECLARE_FUNC(diag_affine_banded_trace, SUFFIX)(
	struct sea_result *aln,
	struct sea_params param,
	char *mat,
	struct mpos o);

sea_int_t
DECLARE_FUNC_GLOBAL(diag_affine_banded, SUFFIX)(
	struct sea_result *aln,
	struct sea_params param,
	void *mat)
{
	sea_int_t retval = SEA_ERROR;
	struct mpos o;
	DECLARE_BENCH(fill);
	DECLARE_BENCH(search);
	DECLARE_BENCH(trace);

	/**
	 * fill in a matrix
	 */
	START_BENCH(fill);
	o = CALL_FUNC(diag_affine_banded_fill, SUFFIX)(
		aln, param,
		(char *)mat + 3 * BYTES_PER_LINE);
	END_BENCH(fill);
	START_BENCH(search);
	o = CALL_FUNC(diag_affine_banded_search, SUFFIX)(
		aln, param,
		(char *)mat + 3 * BYTES_PER_LINE,
		o);
	END_BENCH(search);

#ifdef DEBUG
	fprintf(stderr, "ldmat, %d, %d\n", aln->alen+BAND_WIDTH, aln->blen+BAND_WIDTH);
	for(int i = 0; i < aln->alen+1; i++) {
		for(int j = 0; j < aln->blen+1; j++) {
			fprintf(stderr, "%d, %d, %d, %d, %d, %d\n", i, j, 
            COP(i, j), COQ(i, j), AADDR(COP(i, j), COQ(i, j)), ASCOREV((char *)mat + 
            AADDR(COP(i, j), COQ(i, j))));
		}
	}
	fprintf(stderr, "\n");

	for(int i = 0; i < 3*o.e.p+4; i++) {
		for(int j = 0; j < BAND_WIDTH; j++) {
			fprintf(stderr, "%4d, ", ASCOREV((char *)mat + i*BYTES_PER_LINE + j*BYTES_PER_CELL));
		}
		fprintf(stderr, "\n");
	}
#endif

	
	if(o.m.p == -1) { return SEA_ERROR_OVERFLOW; }
	if(aln->aln == NULL) { return SEA_SUCCESS; }
	START_BENCH(trace);
	retval = CALL_FUNC(diag_affine_banded_trace, SUFFIX)(
		aln, param,
		(char *)mat + 3 * BYTES_PER_LINE,
		o);
	END_BENCH(trace);
	return(retval);
}


struct mpos
DECLARE_FUNC(diag_affine_banded_fill, SUFFIX)(
	struct sea_result *aln,
	struct sea_params param,
	char *mat)
{
	sea_int_t const bw = BAND_WIDTH;
	sea_int_t i, j;
	sea_int_t m = param.m,
			  x = param.x,
			  gi = param.gi,
			  ge = param.ge,
			  xdrop = param.xdrop;
	sea_int_t apos = 0,				/** seq.a read position */
			  bpos = 0;				/** seq.b read position */
	sea_int_t alen = aln->alen,
			  blen = aln->blen;
	sea_int_t alim = alen + bw/2,
			  blim = blen + bw/2;
	struct mpos o = {{0, 0, 0, 0}, {0, 0, 0, 0}};

	DECLARE_VEC_CELL_REG(mv);
	DECLARE_VEC_CELL_REG(xv);
	DECLARE_VEC_CELL_REG(giv);
	DECLARE_VEC_CELL_REG(gev);
	DECLARE_VEC_CHAR_REG(wq);
	DECLARE_VEC_CHAR_REG(wt);
	DECLARE_VEC_CELL_REG(v);
	DECLARE_VEC_CELL_REG(pv);
	DECLARE_VEC_CELL_REG(f);
	DECLARE_VEC_CELL_REG(e);
	DECLARE_VEC_CELL_REG(tmp1);
	DECLARE_VEC_CELL_REG(tmp2);
	DECLARE_VEC_CELL_REG(maxv);
	#if ALG == SW
		DECLARE_VEC_CELL_REG(zv);
	#endif
	DECLARE_SEQ(a);
	DECLARE_SEQ(b);

	CLEAR_SEQ(a, aln->a, aln->apos, aln->alen);
	CLEAR_SEQ(b, aln->b, aln->bpos, aln->blen);

	VEC_SET(mv, m);					
	VEC_SET(xv, x);					
	VEC_SET(giv, gi);				
	VEC_SET(gev, ge);				
	VEC_SET(maxv, CELL_MIN);		
	VEC_SET(pv, CELL_MIN);			
	#if ALG == SW
		VEC_SET(zv, 0);				
	#endif

	VEC_INIT_PVN(v, m, gi, i); VEC_STORE(mat, v);
	VEC_SET(f, CELL_MIN); VEC_STORE(mat, f);
	VEC_SET(e, CELL_MIN); VEC_STORE(mat, e);
	VEC_MAX(maxv, maxv, v);

	VEC_CHAR_SETZERO(wq);			
	VEC_CHAR_SETONES(wt);			
	for(apos = 0; apos < bw/2; apos++) {
		FETCH(a, apos); PUSHQ((apos < alen) ? DECODE(a) : 0,    wq);
	}
	for(bpos = 0; bpos < bw/2-1; bpos++) {
		FETCH(b, bpos); PUSHT((bpos < blen) ? DECODE(b) : 0xff, wt);
	}
	#if ALG == SW
	 	#define VEC_SW_MAX(a, b, c)		VEC_MAX(a, b, c);
	#elif ALG == NW || ALG == SEA || ALG == XSEA
	 	#define VEC_SW_MAX(a, b, c)		;
	#endif
	#if ALG == SW || ALG == SEA || ALG == XSEA
	 	#define VEC_SEA_MAX(a, b, c)	VEC_MAX(a, b, c);
	#elif ALG == NW
	 	#define VEC_SEA_MAX(a, b, c)	;
	#endif
	
	

	i = 0; j = 0;							/** the center cell of the init vector */
	while(i < alim && j < blim) {
		j += 1;								/** go downward */
		UPDATE_FHALF();
		VEC_SHIFT_R(e); VEC_INSERT_MSB(e, CELL_MIN);
		FETCH(b, bpos); PUSHT((bpos < blen) ? DECODE(b) : 0xff, wt); bpos++;
		UPDATE_LHALF();
		VEC_STORE(mat, v); VEC_STORE(mat, f); VEC_STORE(mat, e);

		i += 1;								/** go rightward */
		UPDATE_FHALF();
		VEC_SHIFT_L(f); VEC_INSERT_LSB(f, CELL_MIN);
		FETCH(a, apos); PUSHQ((apos < alen) ? DECODE(a) : 0, wq); apos++;
		UPDATE_LHALF();
		VEC_STORE(mat, v); VEC_STORE(mat, f); VEC_STORE(mat, e);

		if(ALG == XSEA && VEC_CENTER(v) + xdrop - VEC_CENTER(maxv) < 0) { break; }
	}
	#undef VEC_SW_MAX
	#undef VEC_SEA_MAX
	#undef UPDATE_FHALF
	#undef UPDATE_LHALF

	aln->len = COP(alim, blim);
	o.e.i = i; o.e.j = j;
	o.e.p = COP(i, j); o.e.q = COQ(i, j);	/** o.e.q == 0 */
	if(ALG != NW) {
		VEC_STORE(mat, maxv);				/** store max vector at the end of the memory */
		VEC_ASSIGN(tmp1, maxv);
		for(i = 1; i < bw; i++) {
			VEC_SHIFT_R(tmp1);
			VEC_MAX(maxv, tmp1, maxv);		/** extract maximum score in the maxv vector */ 
		}
		VEC_STORE(mat, maxv);				/** store max of the max vector */
	}
	return(o);
}


struct mpos
DECLARE_FUNC(diag_affine_banded_search, SUFFIX)(
	struct sea_result *aln,
	struct sea_params param,
	char *mat,
	struct mpos o)
{
	sea_int_t const bw = BAND_WIDTH,
					bc = BYTES_PER_CELL,
					bl = BYTES_PER_LINE;

	if(ALG == NW) {
		o.m.i = aln->alen; o.m.j = aln->blen;
		o.m.p = COP(aln->alen, aln->blen); o.m.q = COQ(aln->alen, aln->blen);
	} else {	/** ALG == SEA or ALG == SW */
		sea_int_t p, q;
		char *lmat = mat + AADDR(o.e.p+1, -bw/2),
			 *tmat;
		sea_int_t max = ASCOREV(lmat + bl);

		o.m.i = o.m.j = o.m.p = o.m.q = 0;
		if(max == CELL_MAX) { o.m.p = -1; return(o); }
		if(max == 0) { return(o); } /** wind back to (0, 0) if no meaningful score was found */
		for(q = -bw/2; q < bw/2; q++, lmat += bc) {
			if(ASCOREV(lmat) == max) {
				for(p = o.e.p, tmat = lmat-3*bl;
					p >= 0 && ASCOREV(tmat) != max;
					p--, tmat -= 3*bl) {}
				if(p >= o.m.p) { o.m.i = COX(p, q); o.m.j = COY(p, q); o.m.p = p; o.m.q = q; }
			}
		}
		aln->alen = o.m.i;
		aln->blen = o.m.j;
	}
	aln->score = ASCOREV(mat + AADDR(o.m.p, o.m.q));
	return(o);
}


sea_int_t
DECLARE_FUNC(diag_affine_banded_trace, SUFFIX)(
	struct sea_result *aln,
	struct sea_params param,
	char *mat,
	struct mpos o)
{
	sea_int_t mi = o.m.i,
			  mj = o.m.j,
			  mp = o.m.p,
			  mq = o.m.q;
	sea_int_t m = param.m,
			  x = param.x,
			  gi = param.gi,
			  ge = param.ge;
	char *tmat = (char *)mat + AADDR(mp, mq);
	sea_int_t score, cost, diag;

	DECLARE_SEQ(a);
	DECLARE_SEQ(b);
	DECLARE_ALN(l);

	CLEAR_SEQ(a, aln->a, aln->apos, aln->alen);
	CLEAR_SEQ(b, aln->b, aln->bpos, aln->blen);
	CLEAR_ALN(l, aln->aln, aln->len);

	score = ASCOREV(tmat);
	FETCH(a, mi-1); FETCH(b, mj-1);
	while(mi > 0 && mj > 0) {
		diag = ASCOREV(tmat + ATOPLEFT(mp, mq));
		cost = COMPARE(a, b) ? m : x;
		if(score == (diag + cost)) {
			tmat += ATOPLEFT(mp, mq); mp -= 2;
			mi--; FETCH(a, mi-1);
			mj--; FETCH(b, mj-1);
			PUSH(l, (cost == m) ? MATCH_CHAR : MISMATCH_CHAR);
			if(ALG == SW && score <= cost) {
				aln->apos += mi; aln->bpos += mj;
				aln->alen -= mi; aln->blen -= mj;
				mi = mj = 0; break;
			}
			score = diag;
		} else if(score == ASCOREE(tmat)) {
			while(mi > 0 && ASCOREE(tmat) == ASCOREE(tmat + ALEFT(mp, mq)) + ge) {
				tmat += ALEFT(mp, mq); mq += ALEFTQ(mp, mq); mp--;
				mi--; PUSH(l, DELETION_CHAR);
			}
			tmat += ALEFT(mp, mq); mq += ALEFTQ(mp, mq); mp--;
			mi--; FETCH(a, mi-1); PUSH(l, DELETION_CHAR);
			score = ASCOREV(tmat);
		} else if(score == ASCOREF(tmat)) {
			while(mj > 0 && ASCOREF(tmat) == ASCOREF(tmat + ATOP(mp, mq)) + ge) {
				tmat += ATOP(mp, mq); mq += ATOPQ(mp, mq); mp--;
				mj--; PUSH(l, INSERTION_CHAR);
			}
			tmat += ATOP(mp, mq); mq += ATOPQ(mp, mq); mp--;
			mj--; FETCH(b, mj-1); PUSH(l, INSERTION_CHAR);
			score = ASCOREV(tmat);
		} else {
			return SEA_ERROR_OUT_OF_BAND;
		}
	}
	while(mi > 0) { mi--; PUSH(l, DELETION_CHAR); }
	while(mj > 0) { mj--; PUSH(l, INSERTION_CHAR); }
	aln->len = LENGTH(l);
	REVERSE(l);
	return SEA_SUCCESS;
}

sea_int_t
DECLARE_FUNC_GLOBAL(diag_affine_banded_matsize, SUFFIX)(
	sea_int_t alen,
	sea_int_t blen,
	sea_int_t bandwidth)
{
	/** barriar (1) + matrix (alen+blen+bw) + max vectors (2) */
	return(3 * (1 + alen + blen + bandwidth + 2) * BYTES_PER_LINE);
}


#ifdef MAIN

int main(int argc, char *argv[])
{
	int matlen, alnlen;
	void *mat;
	struct sea_result aln;
	struct sea_params param;

	param.flags = 0;
	param.m = 2;
	param.x = -3;
	param.gi = -4;
	param.ge = -1;
	param.xdrop = 12;
	param.bandwidth = 32;

	aln.a = argv[1];
	aln.alen = strlen(aln.a);
	aln.apos = 0;
	aln.b = argv[2];
	aln.blen = strlen(aln.b);
	aln.bpos = 0;
	alnlen = aln.alen + aln.blen;
	aln.len = alnlen;

	alnlen = aln.alen + aln.blen;
	matlen = CALL_FUNC(diag_affine_banded_matsize, SUFFIX)(
		aln.alen, aln.blen, param.bandwidth);

	aln.aln = (void *)malloc(alnlen);
	mat = (void *)malloc(matlen);
	CALL_FUNC(diag_affine_banded, SUFFIX)(&aln, param, mat);

	printf("%d, %d, %s\n", aln.score, aln.len, aln.aln);

	free(mat);
	free(aln.aln);
	return 0;
}

#endif
#ifdef TEST
#if SEQ == ascii && ALN == ascii

extern sea_int_t
DECLARE_FUNC_GLOBAL(naive_affine_banded_matsize, DEFAULT_SUFFIX)(
	sea_int_t alen,
	sea_int_t blen,
	sea_int_t bandwidth);

extern sea_int_t
DECLARE_FUNC_GLOBAL(naive_affine_banded, DEFAULT_SUFFIX)(
	struct sea_result *aln,
	struct sea_params param,
	void *mat);


void
DECLARE_FUNC_GLOBAL(test_8_cross_diag_affine_banded, SUFFIX)(
	void)
{
	int i;
	int const cnt = 5;
	sea_int_t m = 2,
			  x = -3,
			  gi = -4,
			  ge = -1;
	char *a, *b;
	struct sea_context *full, *band;

	#if BIT_WIDTH == 8
		int const len = 50;
	#elif BIT_WIDTH == 16
		int const len = 1000;
	#else
		#error "bit width must be 8 or 16."
	#endif

	full = sea_init_fp(
		SEA_BANDWIDTH_64,
		CALL_FUNC(naive_affine_banded, DEFAULT_SUFFIX),
		CALL_FUNC(naive_affine_banded_matsize, DEFAULT_SUFFIX),
		m, x, gi, ge,
		10000);
	band = sea_init_fp(
		SEA_BANDWIDTH_64,
		CALL_FUNC(diag_affine_banded, SUFFIX),
		CALL_FUNC(diag_affine_banded_matsize, SUFFIX),
		m, x, gi, ge,
		10000);


	for(i = 0; i < cnt; i++) {
		a = rseq(len);
		b = mseq(a, 20, 100, 100);
		sea_assert_cross(
			full, band, a, b);
		free(a); free(b);
	}

	sea_clean(full);
	sea_clean(band);
	return;
}

#endif 

#endif	
	
/