.. index:: cloud; sphinx theme, sphinx theme; cloud

=======================
A simple real business cycle model without a labour-leisure choice
=======================

Introduction
==============

This model was popularized by these people.

Model Setup
============

*Households*

In our simple RBC model involving no labour-leisure choice, our housholds' momentary utility function depends only on consumption in period t.
Also, households are assumed to pick levels of consumption which maximize this function over the their entire (typically infinite) lifetimes.

.. math::
   :label: hh_util_func

   \max_{c_t}U(c_t)=\frac{c_t^{1-\rho}}{1-\rho}

The household's budget contstraint can be expressed in a number of different ways, all of which are equivalent to each other:

.. math::
   :label: hh_budget_con1

   y_t = c_t + i_t

which is the most simple method of stating the houshold's budget constraint, simply implying that output has to be exhausted on consumption
investment.

.. math::
   :label: hh_budget_con2

   e^{z_t}AK^{\alpha}L^{1-\alpha} = c_t + k_{t+1}-\left(1-\delta\right)k_t

The second method more explicitly spells out how output is produced via a CRS Cobb-Douglas production function as well as how investment
is related to the perpetual inventory model.

.. math::
   :label: hh_budget_con3s

   w_tl_t+r_tk_t = c_t+i_t

Finally, the last method describes how total income is exhausted through expenditure on consumption and investment. Since the household optimizes
over her its lifetime subject to the budget constraint, the complete problem would be expressible as:

.. math::
   :label: hh_budget_con3s

   U = \sum_{t=0}^\infty \beta^t\left\{\frac{c_t^{1-\rho}}{1-\rho}-\lambda_t\left(c_t+i_t-y_t\right)\right\}

*Firms*

Firms in this model are assumed to be profit-maximizing and to be renting factors of production in competitive input factor markets. This means
that for one of labour they pay the competitive wage rate :math:`w_t` and for one unit of physical capital they pay the real
interest rate :math:`r_t`. Therefore:

.. math::
   :label: firm_profit

   \max_{l_t,k_t}\Pi\left(l_t,k_t\right) = e^{z_t}Ak_t^{1-\alpha}l_t^{\alpha}-w_tl_t-r_tk_t

Strictly speaking we could also spell out the firm's problem as one which is solved over infinitely many periods. However, since the firm
faces an identical problem in each time period, there is no inter-temporality involved here, so we can just focus on the within-period problem
which would be optimal for all periods.

First-Order Conditions of Optimality
====================================

The first-order conditions of optimality are simply obtained by setting up both the household's and the firm's contrained optimisation problems
and taking first derivatives.

*Households*

.. math::
   :label: hh_foc_c

   c_t:\quad\frac{\partial U(c_t)}{\partial c_t}-\lambda_t = 0 \quad \Longleftrightarrow \quad c_t^{-\rho}=\lambda_t

where :math:`\lambda_t` is simply equal to the shadow price of one unit of wealth. This condition simply states that at an optimum
the marginal utility of consumption has to equal the marginal value in utility terms of one extra unit of wealth. Households also have to decide
on how much of their wealth to invest in the physical capital storage technology :math:`k_t`, formally a decision of how much to save:

.. math::
   :label: hh_foc_k

   k_t:\quad\beta\lambda_{t+1}(1+r_{t+1})-\lambda_t = 0

When using the FOC for consumption as well as being more explicit about the real rate of return, we can also write the above as [#f1]_ :

.. math::
   :label: hh_foc_k2

   k_{t+1}:\quad\beta\frac{\partial U(c_{t+1})}{\partial c_{t+1}}(1+\frac{\partial F(l_{t+1},k_{t+1})}{\partial k_{t+1}}-\delta)-\frac{\partial U(c_t)}{\partial c_t} = 0

*Firms*

Firms have to choose optimal quantities of labour and physical capital in order to produce their output and maximize their profits. This leads
to the first-order conditions of optimality:

.. math::
   :label: firm_foc_lab

   l_t:\quad\frac{\partial F(k_t,l_t)}{\partial l_t} - w_t = 0 \Longleftrightarrow e^{z_t}Ak_t^{1-\alpha}l_t^{\alpha-1} = w_t

and for physical capital:

.. math::
   :label: firm_foc_lab

   k_t:\quad\frac{\partial F(k_t,l_t)}{\partial k_t} - r_t = 0 \Longleftrightarrow e^{z_t}Ak_t^{-\alpha}l_t^{\alpha} = r_t

any solution needs to respect the original budget constraint:

.. math::
   :label: lambda_foc

   \lambda_t:\quad y_t-c_t-i_t = 0 \Longleftrightarrow y_t = c_t + i_t


.. rubric:: Footnotes

.. [#f1] Journal articles and text book treatments often use different notations for next-period capital. Sometimes it is written as
         :math:`k_t` to stress the fact that next-period :math:`t+1` capital is determined in this period :math:`t`, while at other times
         it is written as :math:`k_{t+1}` to stress the fact that this will be the amount of capital available next period after it
         was determined in this period.