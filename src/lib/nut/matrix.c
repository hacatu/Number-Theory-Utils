#include "nut/modular_math.h"
#include <nut/matrix.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

bool nut_i64_Matrix_init(nut_i64_Matrix *self, int64_t rows, int64_t cols){
	size_t nbytes;
	if(rows < 0 || cols < 0 || __builtin_mul_overflow(rows, cols, &nbytes) || __builtin_mul_overflow(nbytes, sizeof(int64_t), &nbytes)){
		return false;
	}
	if(!(self->buf = malloc(nbytes))){
		return false;
	}
	self->rows = rows;
	self->cols = cols;
	return true;
}

void nut_i64_Matrix_copy(nut_i64_Matrix *restrict dest, const nut_i64_Matrix *restrict src){
	memcpy(dest->buf, src->buf, src->rows*src->cols*sizeof(int64_t));
}

void nut_i64_Matrix_destroy(nut_i64_Matrix *self){
	free(self->buf);
	self->buf = NULL;
}

void nut_i64_Matrix_fprint(const nut_i64_Matrix *restrict self, FILE *file){
	fprintf(file, "[");
	for(int64_t i = 0; i < self->rows; ++i){
		fprintf(file, "[");
		for(int64_t j = 0; j < self->cols; ++j){
			int64_t a = self->buf[self->cols*i + j];
			if(j + 1 == self->cols){
				fprintf(file, "%"PRIi64, a);
			}else{
				fprintf(file, "%"PRIi64", ", a);
			}
		}
		fprintf(file, "]");
		if(i + 1 < self->cols){
			fprintf(file, ",\n");
		}
	}
	fprintf(file, "]");
}

void nut_i64_Matrix_mul_vec(const nut_i64_Matrix *restrict self, const int64_t vec[restrict static self->cols], int64_t out[restrict static self->rows]){
	if(self->rows){
		memset(out, 0, self->rows*sizeof(int64_t));
	}
	for(int64_t row = 0; row < self->rows; ++row){
		for(int64_t col = 0; col < self->cols; ++col){
			out[row] += self->buf[self->cols*row + col]*vec[col];
		}
	}
}

bool nut_i64_Matrix_fill_I(nut_i64_Matrix *self){
	if(self->rows != self->cols){
		return false;
	}
	if(!self->rows){
		return true;
	}
	memset(self->buf, 0, self->rows*self->cols*sizeof(int64_t));
	for(int64_t i = 0; i < self->rows; ++i){
		self->buf[self->rows*i + i] = 1;
	}
	return true;
}

void nut_i64_Matrix_fill_short_pascal(nut_i64_Matrix *self){
	if(self->rows*self->cols == 0){
		return;
	}
	for(int64_t row = 0; row < self->rows; ++row){
		self->buf[self->cols*row] = 1;
		if(row){
			self->buf[self->cols*row + 1] = row + 1;
		}
		for(int64_t col = 2; col < row && col < self->cols; ++col){
			self->buf[self->cols*row + col] = self->buf[self->cols*(row - 1) + col - 1] + self->buf[self->cols*(row - 1) + col];
		}
		if(row < self->cols){
			self->buf[self->cols*row + row] = row + 1;
		}
		if(row + 1 < self->cols){
			memset(self->buf + self->cols*row + row + 1, 0, (self->cols - row - 1)*sizeof(int64_t));
		}
	}
}

bool nut_i64_Matrix_scale_row(nut_i64_Matrix *self, int64_t row, int64_t col_start, int64_t a){
	for(int64_t i = col_start; i < self->cols; ++i){
		self->buf[self->cols*row + i] *= a;
	}
	return true;
}

bool nut_i64_Matrix_addmul_row(nut_i64_Matrix *self, int64_t i, int64_t j, int64_t a){
	for(int64_t k = 0; k < self->cols; ++k){
		self->buf[self->cols*j + k] += a*self->buf[self->cols*i + k];
	}
	return true;
}

int64_t nut_i64_Matrix_invert_ltr(nut_i64_Matrix *restrict self, nut_i64_Matrix *restrict out){
	if(self->rows != self->cols || self->rows != out->rows || self->cols != out->rows){
		return 0;
	}
	for(int64_t i = 0; i < self->rows; ++i){
		if(!self->buf[self->cols*i + i]){
			return 0;
		}
	}
	int64_t denom = 1;
	nut_i64_Matrix_fill_I(out);
	for(int64_t i = 0; i < self->rows; ++i){
		int64_t a = self->buf[self->cols*i + i];
		for(int64_t j = i + 1; j < self->rows; ++j){
			int64_t b = self->buf[self->cols*j + i];
			int64_t g = nut_i64_egcd(a, b, NULL, NULL);
			if(a != g){
				nut_i64_Matrix_scale_row(self, j, i, a/g);
				nut_i64_Matrix_scale_row(out, j, 0, a/g);
			}
			if(b){
				nut_i64_Matrix_addmul_row(out, i, j, -b/g);
			}
		}
	}
	for(int64_t i = 0; i < self->rows; ++i){
		denom = nut_i64_lcm(denom, self->buf[self->cols*i + i]);
	}
	for(int64_t i = 0; i < self->rows; ++i){
		nut_i64_Matrix_scale_row(out, i, 0, denom/self->buf[self->cols*i + i]);
	}
	return denom;
}

void nut_i64_Matrix_fill_vandemond_vec(uint64_t x, uint64_t k, uint64_t m, int64_t out[static k + 1]){
	if(m){
		x %= m;
	}
	for(uint64_t e = 1, xe = x; e <= k; ++e, xe = m ? xe*x%m : xe*x){
		out[e - 1] = xe;
	}
}

bool nut_u64_ModMatrix_init(nut_u64_ModMatrix *self, uint64_t rows, uint64_t cols, uint64_t modulus){
	size_t nbytes;
	if(__builtin_mul_overflow(rows, cols, &nbytes) || __builtin_mul_overflow(nbytes, sizeof(uint64_t), &nbytes)){
		return false;
	}
	if(!(self->buf = malloc(nbytes))){
		return false;
	}
	self->rows = rows;
	self->cols = cols;
	self->modulus = modulus;
	return true;
}

void nut_u64_ModMatrix_copy(nut_u64_ModMatrix *restrict dest, const nut_u64_ModMatrix *restrict src){
	memcpy(dest->buf, src->buf, src->rows*src->cols*sizeof(uint64_t));
}

void nut_u64_ModMatrix_destroy(nut_u64_ModMatrix *self){
	free(self->buf);
	self->buf = NULL;
}

void nut_u64_ModMatrix_fprint(const nut_u64_ModMatrix *restrict self, FILE *file){
	fprintf(file, "[");
	for(uint64_t i = 0; i < self->rows; ++i){
		fprintf(file, "[");
		for(uint64_t j = 0; j < self->cols; ++j){
			uint64_t a = self->buf[self->cols*i + j];
			if(j + 1 == self->cols){
				fprintf(file, "%"PRIu64, a);
			}else{
				fprintf(file, "%"PRIu64", ", a);
			}
		}
		fprintf(file, "]");
		if(i + 1 < self->cols){
			fprintf(file, ",\n");
		}
	}
	fprintf(file, "] mod %"PRIu64"", self->modulus);
}

bool nut_u64_ModMatrix_fill_I(nut_u64_ModMatrix *self){
	if(self->rows != self->cols){
		return false;
	}
	if(!self->rows){
		return true;
	}
	memset(self->buf, 0, self->rows*self->cols*sizeof(uint64_t));
	for(uint64_t i = 0; i < self->rows; ++i){
		self->buf[self->rows*i + i] = 1;
	}
	return true;
}

void nut_u64_ModMatrix_mul_vec(const nut_u64_ModMatrix *restrict self, const uint64_t vec[restrict static self->cols], uint64_t out[restrict static self->rows]){
	if(self->rows){
		memset(out, 0, self->rows*sizeof(uint64_t));
	}
	for(uint64_t row = 0; row < self->rows; ++row){
		for(uint64_t col = 0; col < self->cols; ++col){
			out[row] = (out[row] + self->buf[self->cols*row + col]*vec[col])%self->modulus;
		}
	}
}

void nut_u64_ModMatrix_fill_short_pascal(nut_u64_ModMatrix *self){
	if(self->rows*self->cols == 0){
		return;
	}
	for(uint64_t row = 0; row < self->rows; ++row){
		self->buf[self->cols*row] = 1;
		if(row){
			self->buf[self->cols*row + 1] = (row + 1)%self->modulus;
		}
		for(uint64_t col = 2; col < row && col < self->cols; ++col){
			uint64_t v = self->buf[self->cols*(row - 1) + col - 1] + self->buf[self->cols*(row - 1) + col];
			self->buf[self->cols*row + col] = v >= self->modulus ? v - self->modulus : v;
		}
		if(row < self->cols){
			self->buf[self->cols*row + row] = (row + 1)%self->modulus;
		}
		if(row + 1 < self->cols){
			memset(self->buf + self->cols*row + row + 1, 0, (self->cols - row - 1)*sizeof(uint64_t));
		}
	}
}

bool nut_u64_ModMatrix_scale_row(nut_u64_ModMatrix *self, uint64_t row, uint64_t col_start, uint64_t a){
	for(uint64_t i = col_start; i < self->cols; ++i){
		self->buf[self->cols*row + i] *= a;
		self->buf[self->cols*row + i] %= self->modulus;
	}
	return true;
}

bool nut_u64_ModMatrix_addmul_row(nut_u64_ModMatrix *self, uint64_t i, uint64_t j, uint64_t a){
	for(uint64_t k = 0; k < self->cols; ++k){
		self->buf[self->cols*j + k] += a*self->buf[self->cols*i + k];
		self->buf[self->cols*j + k] %= self->modulus;
	}
	return true;
}

bool nut_u64_ModMatrix_invert_ltr(nut_u64_ModMatrix *restrict self, nut_u64_ModMatrix *restrict out){
	if(self->rows != self->cols || self->rows != out->rows || self->cols != out->rows){
		return 0;
	}
	for(uint64_t i = 0; i < self->rows; ++i){
		if(!self->buf[self->cols*i + i]){
			return false;
		}
	}
	nut_u64_ModMatrix_fill_I(out);
	for(uint64_t i = 0; i < self->rows; ++i){
		uint64_t a = self->buf[self->cols*i + i];
		int64_t _ainv = nut_i64_modinv(a, self->modulus);
		nut_u64_ModMatrix_scale_row(out, i, 0, (uint64_t)_ainv);
		for(uint64_t j = i + 1; j < self->rows; ++j){
			uint64_t b = self->buf[self->cols*j + i];
			if(b){
				nut_u64_ModMatrix_addmul_row(out, i, j, self->modulus - b);
			}
		}
	}
	return true;
}

void nut_u64_Matrix_fill_vandemond_vec(uint64_t x, uint64_t k, uint64_t m, uint64_t out[static k + 1]){
	for(uint64_t e = 1, xe = x; e <= k; ++e, xe = xe*x%m){
		out[e - 1] = xe;
	}
}

