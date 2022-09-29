	# Compute (a + b + carry_in) and store the result in the address
	# of sum, the value 0 is stored at the address of carry_out unless
        # the expression above produces a carry in which case the
	# value 1 is written.
	# %rdi: a
        # %rsi: b
        # %rdx: carry_in
        # %rcx: address of sum
        # %r8:  address of carry_out
        .globl add
add:
        # Add your code here ...
	movq %rdx, (%rcx)
	addq %rdi, (%rcx)
	jnc next
	movq $1, (%r8)
	jmp next2
next2:	addq %rsi, (%rcx)
	jmp end
next:	addq %rsi, (%rcx)
	jnc nocarr
	movq $1, (%r8)
	jmp end
nocarr: movq $0, (%r8)
end:	retq
