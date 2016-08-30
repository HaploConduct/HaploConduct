/* static_bitsequence_tester.h
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * static_bitsequence_tester
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <sys/stat.h>
#include <iostream>
#include <sstream>
#include <basics.h>
#include <static_bitsequence.h>
#include <alphabet_mapper.h>
#include <static_sequence.h>
#include <static_sequence_builder.h>


#ifndef STATIC_BITSEQUENCE_TESTER_H
#define STATIC_BITSEQUENCE_TESTER_H

void load(char *fname, uint ** text, uint * n);
void test_bitsequence(uint * bitseq, uint len, static_bitsequence * bs);
void speed_access(static_bitsequence * ss, uint * bitseq, uint n);
void speed_rank0(static_bitsequence * ss, uint * bitseq, uint n);
void speed_rank1(static_bitsequence * ss, uint * bitseq, uint n);
void speed_select0(static_bitsequence * ss, uint * bitseq, uint n);
void speed_select1(static_bitsequence * ss, uint * bitseq, uint n);
void speed_selectnext1(static_bitsequence * ss, uint * bitseq, uint n);

#endif
