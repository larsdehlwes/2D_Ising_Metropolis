template <typename Container, typename T = typename std::decay<decltype(*std::begin(std::declval<Container>()))>::type> T stdev(Container && c)
{ 
    auto b = std::begin(c), e = std::end(c);
    auto size = std::distance(b, e);
    auto sum = std::accumulate(b, e, T());
    auto mean = sum / size;
    T accum = T();
    for (const auto d : c)
        accum += (d - mean) * (d - mean);
    return std::sqrt(accum / (size - 1));
}

template <typename Container, typename T = typename std::decay<decltype(*std::begin(std::declval<Container>()))>::type> T avg(Container && c)
{
    auto b = std::begin(c), e = std::end(c);
    auto size = std::distance(b, e);
    auto sum = std::accumulate(b, e, T());
    auto mean = sum / size;
    return mean;
}

double corr(const std::vector<double>& x, const std::vector<double>& y) {
    const auto n    = x.size();
    const auto s_x  = std::accumulate(x.begin(), x.end(), 0.0);
    const auto s_y  = std::accumulate(y.begin(), y.end(), 0.0);
    const auto s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
    const auto s_yy = std::inner_product(y.begin(), y.end(), y.begin(), 0.0);
    const auto s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
    const auto a    = (n * s_xy - s_x * s_y) / ((n * s_xx - s_x * s_x) * (n * s_yy - s_y * s_y));
    return a;
}

double slope(const std::vector<double>& x, const std::vector<double>& y) {
    const auto n    = x.size();
    const auto s_x  = std::accumulate(x.begin(), x.end(), 0.0);
    const auto s_y  = std::accumulate(y.begin(), y.end(), 0.0);
    const auto s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
    const auto s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
    const auto a    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
    return a;
}
