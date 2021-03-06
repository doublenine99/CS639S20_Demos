# mark_description "Intel(R) C Intel(R) 64 Compiler for applications running on Intel(R) 64, Version 19.1.0.166 Build 20191121";
# mark_description "-S -xCOMMON-AVX512 -fno-alias";
	.file "MatMatMultiplyBlockHelper.cpp"
	.text
..TXTST0:
.L_2__routine_start__Z25MatMatMultiplyBlockHelperRA64_A64_KfS2_RA64_A64_f_0:
# -- Begin  _Z25MatMatMultiplyBlockHelperRA64_A64_KfS2_RA64_A64_f
	.text
# mark_begin;
       .align    16,0x90
	.globl _Z25MatMatMultiplyBlockHelperRA64_A64_KfS2_RA64_A64_f
# --- MatMatMultiplyBlockHelper(const float (&)[64][64], const float (&)[64][64], float (&)[64][64])
_Z25MatMatMultiplyBlockHelperRA64_A64_KfS2_RA64_A64_f:
# parameter 1: %rdi
# parameter 2: %rsi
# parameter 3: %rdx
..B1.1:                         # Preds ..B1.0
                                # Execution count [1.00e+00]
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
..___tag_value__Z25MatMatMultiplyBlockHelperRA64_A64_KfS2_RA64_A64_f.1:
..L2:
                                                          #6.1
        xorb      %cl, %cl                                      #9.5
        xorl      %eax, %eax                                    #9.5
                                # LOE rax rdx rbx rbp rsi rdi r12 r13 r14 r15 cl
..B1.2:                         # Preds ..B1.4 ..B1.1
                                # Execution count [6.40e+01]
        xorl      %r10d, %r10d                                  #10.5
        lea       (%rdi,%rax), %r9                              #13.40
        vmovups   192(%rax,%rdx), %zmm3                         #14.41
        xorl      %r8d, %r8d                                    #10.5
        vmovups   128(%rax,%rdx), %zmm2                         #14.41
        vmovups   64(%rax,%rdx), %zmm1                          #14.41
        vmovups   (%rax,%rdx), %zmm0                            #14.41
                                # LOE rax rdx rbx rbp rsi rdi r8 r9 r10 r12 r13 r14 r15 cl zmm0 zmm1 zmm2 zmm3
..B1.3:                         # Preds ..B1.3 ..B1.2
                                # Execution count [4.10e+03]
        vbroadcastss (%r9,%r10,4), %zmm4                        #13.40
        incq      %r10                                          #10.5
        vfmadd231ps (%r8,%rsi), %zmm4, %zmm0                    #15.18
        vfmadd231ps 64(%r8,%rsi), %zmm4, %zmm1                  #15.18
        vfmadd231ps 128(%r8,%rsi), %zmm4, %zmm2                 #15.18
        vfmadd231ps 192(%r8,%rsi), %zmm4, %zmm3                 #15.18
        addq      $256, %r8                                     #10.5
        cmpq      $64, %r10                                     #10.5
        jb        ..B1.3        # Prob 98%                      #10.5
                                # LOE rax rdx rbx rbp rsi rdi r8 r9 r10 r12 r13 r14 r15 cl zmm0 zmm1 zmm2 zmm3
..B1.4:                         # Preds ..B1.3
                                # Execution count [6.40e+01]
        incb      %cl                                           #9.5
        vmovups   %zmm3, 192(%rax,%rdx)                         #16.30
        vmovups   %zmm2, 128(%rax,%rdx)                         #16.30
        vmovups   %zmm1, 64(%rax,%rdx)                          #16.30
        vmovups   %zmm0, (%rax,%rdx)                            #16.30
        addq      $256, %rax                                    #9.5
        cmpb      $64, %cl                                      #9.5
        jb        ..B1.2        # Prob 98%                      #9.5
                                # LOE rax rdx rbx rbp rsi rdi r12 r13 r14 r15 cl
..B1.5:                         # Preds ..B1.4
                                # Execution count [1.00e+00]
        vzeroupper                                              #18.1
        ret                                                     #18.1
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	_Z25MatMatMultiplyBlockHelperRA64_A64_KfS2_RA64_A64_f,@function
	.size	_Z25MatMatMultiplyBlockHelperRA64_A64_KfS2_RA64_A64_f,.-_Z25MatMatMultiplyBlockHelperRA64_A64_KfS2_RA64_A64_f
..LN_Z25MatMatMultiplyBlockHelperRA64_A64_KfS2_RA64_A64_f.0:
	.data
# -- End  _Z25MatMatMultiplyBlockHelperRA64_A64_KfS2_RA64_A64_f
	.data
	.section .note.GNU-stack, ""
# End
