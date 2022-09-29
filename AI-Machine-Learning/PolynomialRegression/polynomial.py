import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from statistics import mean

# @author Skye Russ


def generate_coefficients():
    """
    Returns a sequence of five numbers, to be used as coefficients of a polynomial. Each number is chosen uniformly from the
    interval [-0.5, 0.5).
    """
    return [np.random.uniform(-0.5, 0.5) for i in range(5)]


def generate_data(m, coefficients):
    """
    Returns two arrays, X and y, each of which is m by 1. The values of X are evenly spaced across the interval
    [-5.0, 5.0]. For each x, the corresponding y value is

    a + b * X + c * X**2 + d * X**3 + e * X**4 + <noise>

    where coefficients is (a, b, c, d, e) and the noise for each point is normally distributed with mean 0 and
    standard deviation 1.
    """
    x = np.linspace(-5.0, 5.0, m).reshape(m, 1)
    y = np.array(coefficients[0] + coefficients[1] * x + coefficients[2] * (x ** 2) + coefficients[3] * (x ** 3) +
                 coefficients[4] * (x ** 4) + np.random.standard_normal())
    return x, y


def plot_data(X, y):
    """
    Plots X and y (as a scatter plot) and also constrains the y limit so that later, much larger values of y will not
    reset it.
    """
    plt.ylim((y.min() - 0.1 * (y.max() - y.min()),
              y.max() + 0.1 * (y.max() - y.min())))
    plt.scatter(X, y)


def generate_and_plot_data(m):
    """
    Generates m data points and plots them.
    """
    c = generate_coefficients()
    X, y = generate_data(m, c)
    plot_data(X, y)
    plt.show()


def fit_curve(X, y, degree):
    """
    Returns a trained model that fits a polynomial of the specified degree to the data.
    """
    poly_features = PolynomialFeatures(degree)
    X_poly = poly_features.fit_transform(X)
    lin_reg = LinearRegression()
    lin_reg.fit(X_poly, y)
    return lin_reg


def plot_curve(degree, model):
    """
    Plots a curve for model, which represents a polynomial of the specified degree.
    The x values for the curve are 100 points evenly spaced across the interval [-5.0, 5.0].
    """
    # initialize x with uniformly distributed and y with intercept value
    x = np.linspace(-5.0, 5.0, 100)
    y = [model.intercept_[0] for i in range(100)]
    # [0,limit] (range naturally goes from [0,limit)) so we must add 1 to degree
    for i in range(degree + 1):
        y += (x ** i) * model.coef_[0][i]
    plt.plot(x, y, label=f"degree {degree}")


def mse(X, y, degree, model):
    """
    Returns the mean squared error for model (a polynomial of the specified degree) on X and y.
    """
    # initialize y with the intercept
    y_pred = [model.intercept_[0] for i in range(len(X))]
    # [0,limit] (range naturally goes from [0,limit)) so we must add 1 to degree
    for i in range(degree+1):
        for j in range(len(X)):
            y_pred[j] += (X[j][0] ** i) * model.coef_[0][i]
    return mean_squared_error(y, y_pred)


def experiment_1(m):
    """
    Generates m training points and fits models of degrees 1, 2, and 20. Plots the data and the curves for the models.
    """
    coeffs = generate_coefficients()
    X, y = generate_data(m, coeffs)
    plot_data(X, y)
    for d in [1, 2, 7,8,9,20]:
        model = fit_curve(X, y, d)
        plot_curve(d, model)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()


def experiment_2(m):
    """
    Runs the following experiment 100 times:

    Generate m training data points
    Generate 100 testing data points (using the same coefficients)
    For each d from 1 through 30, fit a curve of degree d to the training data and measure its mse on the testing data.

    After the 100 runs, plots the average mse of each degree.
    """
    mses = {i: [] for i in range(1, 31)}
    for i in range(100):
        print(f"Testing {i + 1} of 100")
        coeffs = generate_coefficients()
        X_train, y_train = generate_data(m, coeffs)
        X_test, y_test = generate_data(100, coeffs)
        for d in range(1, 31):
            model = fit_curve(X_train, y_train, d)
            mses[d] += [mse(X_test, y_test, d, model)]
    averages = [mean(mses[d]) for d in mses]
    plt.ylim(0, 500)
    plt.plot(range(1, 31), averages)
    plt.xlabel('Degree')
    plt.ylabel('Average MSE (100 runs)')
    plt.show()


if __name__ == '__main__':
    # generate_and_plot_data(100)
    experiment_1(5)
    experiment_1(20)
    experiment_1(50)
    experiment_1(100)
    experiment_1(200)
    # experiment_2(20)
    pass
