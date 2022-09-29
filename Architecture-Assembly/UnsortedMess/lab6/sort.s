	.file	"sort.c"
	.text
	.globl	sort
	.type	sort, @function
sort:
.LFB0:
	.cfi_startproc
	subq	$1, %rsi
	jmp	.L2
.L3:
	addq	$1, %rax
.L5:
	cmpq	%rax, %rsi
	jle	.L7
	leaq	8(%rdi,%rax,8), %r8
	movq	(%r8), %rcx
	leaq	(%rdi,%rax,8), %rdx
	movq	(%rdx), %r9
	cmpq	%r9, %rcx
	jge	.L3
	movq	%r9, (%r8)
	movq	%rcx, (%rdx)
	jmp	.L3
.L7:
	subq	$1, %rsi
.L2:
	testq	%rsi, %rsi
	jle	.L8
	movl	$0, %eax
	jmp	.L5
.L8:
	rep ret
	.cfi_endproc
.LFE0:
	.size	sort, .-sort
	.ident	"GCC: (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0"
	.section	.note.GNU-stack,"",@progbits
