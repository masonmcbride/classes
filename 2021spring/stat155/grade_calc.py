hw = 0.15666666666666668
mt1 = 0.23
def g(final):
    return max((.80 * final/100) + hw, (.55 * final/100) + mt1 + hw)

def grade_for_a(grade):
    return min([final for final in range(0, 101) if g(final) > (grade/100)])

if __name__ == '__main__':
    print("what grade do you need?")
    grade = int(input("> "))
    print(f"grade needed: {grade_for_a(grade)}")
