double scalar_product_simple(double vector_a[], double vector_b[])
{
    double answer = 0;
    for (int i = 0; i < 3; ++i)
        answer = answer + vector_a[i]*vector_b[i];
    return answer;
}
