import numpy as np
from pytest import approx, raises


def test_make_2x2_array():
    from solcore.absorption_calculator.tmm_core_vec import make_2x2_array
    z = np.zeros(3)
    o = np.ones(3)
    result = make_2x2_array(z, o, z, o, dtype=float)

    assert result == approx(np.array([[[0., 1.], [0., 1.]],
                              [[0., 1.], [0., 1.]],
                              [[0., 1.], [0., 1.]]]))


def test_snell():
    from solcore.absorption_calculator.tmm_core_vec import snell
    assert snell(3, 2, 0.7) == approx(1.3105496419558818)


def test_list_snell():
    from solcore.absorption_calculator.tmm_core_vec import list_snell
    assert list_snell(np.array([2, 3, 4]), 0.3) == approx(np.array([0.3, 0.19831075, 0.14830313]))


def test_interface_r():
    from solcore.absorption_calculator.tmm_core_vec import interface_r
    th1 = 0.2
    n1 = 2
    n2 = 3
    th2 = np.arcsin((n1/n2)*np.sin(th1))

    s_res = interface_r('s', n1, n2, th1, th2)
    p_res = interface_r('p', n1, n2, th1, th2)

    assert [s_res, p_res] == approx([-0.20541108217641596, 0.19457669033430525])

    with raises(ValueError):
        interface_r('not_s_or_p', n1, n2, th1, th2)


def test_interface_t():
    from solcore.absorption_calculator.tmm_core_vec import interface_t
    th1 = 0.2
    n1 = 2
    n2 = 3
    th2 = np.arcsin((n1/n2)*np.sin(th1))

    s_res = interface_t('s', n1, n2, th1, th2)
    p_res = interface_t('p', n1, n2, th1, th2)

    assert [s_res, p_res] == approx([0.7945889178235841, 0.7963844602228701])

    with raises(ValueError):
        interface_t('not_s_or_p', n1, n2, th1, th2)


def test_R_from_r():
    from solcore.absorption_calculator.tmm_core_vec import R_from_r
    assert R_from_r(np.sqrt(2) + 1j*np.sqrt(2)) == 4.0


def test_T_from_t():
    from solcore.absorption_calculator.tmm_core_vec import T_from_t
    th1 = 0.2
    n1 = 2
    n2 = 3
    th2 = np.arcsin((n1/n2)*np.sin(th1))
    ts = 0.7945889178235841
    tp = 0.7963844602228701
    assert T_from_t('s', ts, n1, n2, th1, th2) == approx(0.9578062873191139)
    assert T_from_t('p', tp, n1, n2, th1, th2) == approx(0.962139911578548)
    with raises(ValueError):
        T_from_t('not_s_or_p', ts, n1, n2, th1, th2)


def test_power_entering_from_r():
    from solcore.absorption_calculator.tmm_core_vec import power_entering_from_r
    rs = -0.20541108217641596
    rp = 0.19457669033430525

    n1 = 2
    th1 = 0.2

    assert power_entering_from_r('s', rs, n1, th1) == approx(0.9578062873191138)
    assert power_entering_from_r('p', rp, n1, th1) == approx(0.962139911578548)

    with raises(ValueError):
        power_entering_from_r('not_s_or_p', rs, n1, th1)


def test_interface_R():
    from solcore.absorption_calculator.tmm_core_vec import interface_R
    th1 = 0.2
    n1 = 2
    n2 = 3
    th2 = np.arcsin((n1/n2)*np.sin(th1))

    assert(interface_R('s', n1, n2, th1, th2) == approx(0.042193712680886314))
    assert(interface_R('p', n1, n2, th1, th2) == approx(0.037860088421452116))


def test_interface_T():
    from solcore.absorption_calculator.tmm_core_vec import interface_T
    th1 = 0.2
    n1 = 2
    n2 = 3
    th2 = np.arcsin((n1/n2)*np.sin(th1))

    assert(interface_T('s', n1, n2, th1, th2) == approx(0.9578062873191139))
    assert(interface_T('p', n1, n2, th1, th2) == approx(0.962139911578548))

### Test for coh_tmm
def test_coh_tmm_exceptions():
    from solcore.absorption_calculator.tmm_core_vec import coh_tmm
    n_list = np.array([[1.5, 1 + 0.4*1j, 2 + 3*1j, 5, 4+1*1j],
                      [1.3, 1.2 + 0.2*1j, 1.5 + 0.3*1j, 4, 3+0.1*1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    with raises(ValueError):
        coh_tmm('s', n_list, d_list, np.array([0.3, 0, 0.8]), lam_vac)

    with raises(ValueError):
        coh_tmm('s', n_list[0], d_list, th_0, lam_vac)

    with raises(ValueError):
        coh_tmm('s', n_list, np.array([10, 200, 187.3, 1973.5, np.inf]), th_0, lam_vac)

    with raises(ValueError):
        coh_tmm('s', n_list, d_list, 0.2+0.2*1j, lam_vac)


def test_coh_tmm_s_r():
    from solcore.absorption_calculator.tmm_core_vec import coh_tmm
    n_list = np.array([[1.5, 1 + 0.4*1j, 2 + 3*1j, 5, 4+1*1j],
                      [1.3, 1.2 + 0.2*1j, 1.5 + 0.3*1j, 4, 3+0.1*1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    result = coh_tmm('s', n_list, d_list, th_0, lam_vac)

    assert result['r'] == approx(np.array([0.14017645-0.2132843j , 0.22307786-0.10704008j]))


def test_coh_tmm_s_t():
    from solcore.absorption_calculator.tmm_core_vec import coh_tmm
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    result = coh_tmm('s', n_list, d_list, th_0, lam_vac)

    assert result['t'] == approx(np.array([ 1.78669633e-05-9.79824244e-06j, -8.86075993e-02-4.05953564e-01j]))


def test_coh_tmm_s_R():
    from solcore.absorption_calculator.tmm_core_vec import coh_tmm
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    result = coh_tmm('s', n_list, d_list, th_0, lam_vac)

    assert result['R'] == approx(np.array([0.06513963, 0.06122131]))


def test_coh_tmm_s_T():
    from solcore.absorption_calculator.tmm_core_vec import coh_tmm
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    result = coh_tmm('s', n_list, d_list, th_0, lam_vac)

    assert result['T'] == approx(np.array([1.15234466e-09, 4.13619185e-01]))


def test_coh_tmm_s_power_entering():
    from solcore.absorption_calculator.tmm_core_vec import coh_tmm
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    result = coh_tmm('s', n_list, d_list, th_0, lam_vac)

    assert result['power_entering'] == approx(np.array([0.93486037, 0.93877869]))


def test_coh_tmm_s_vw_list():
    from solcore.absorption_calculator.tmm_core_vec import coh_tmm
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    result = coh_tmm('s', n_list, d_list, th_0, lam_vac)

    assert result['vw_list'] == approx(np.array([[[0.00000000e+00+0.00000000e+00j,
        0.00000000e+00+0.00000000e+00j],
        [0.00000000e+00+0.00000000e+00j,
        0.00000000e+00+0.00000000e+00j]],
       [[1.18358724e+00-2.33272105e-01j,
         -4.34107939e-02+1.99878010e-02j],
        [1.03160316e+00-7.28921467e-02j,
        1.91474694e-01-3.41479380e-02j]],
       [[-8.59535500e-02+1.06568462e-01j,
         -1.36521327e-09+2.83859953e-10j],
        [6.08369346e-01+5.06683493e-01j,
        1.75320349e-01-9.58306162e-02j]],
       [[-1.23112929e-05+1.37276841e-05j,
         -1.94390395e-06+2.16097082e-06j],
        [-6.54156818e-02+3.57104644e-01j,
        3.38453387e-02+4.04808706e-02j]],
       [[1.78669633e-05-9.79824244e-06j,
        0.00000000e+00+0.00000000e+00j],
        [-8.86075993e-02-4.05953564e-01j,
        0.00000000e+00+0.00000000e+00j]]]))


def test_coh_tmm_p_r():
    from solcore.absorption_calculator.tmm_core_vec import coh_tmm
    n_list = np.array([[1.5, 1 + 0.4*1j, 2 + 3*1j, 5, 4+1*1j],
                      [1.3, 1.2 + 0.2*1j, 1.5 + 0.3*1j, 4, 3+0.1*1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    result = coh_tmm('s', n_list, d_list, th_0, lam_vac)

    assert result['r'] == approx(np.array([0.14017645-0.2132843j , 0.22307786-0.10704008j]))


def test_coh_tmm_p_t():
    from solcore.absorption_calculator.tmm_core_vec import coh_tmm
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    result = coh_tmm('s', n_list, d_list, th_0, lam_vac)

    assert result['t'] == approx(np.array([ 1.78669633e-05-9.79824244e-06j, -8.86075993e-02-4.05953564e-01j]))


def test_coh_tmm_p_R():
    from solcore.absorption_calculator.tmm_core_vec import coh_tmm
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    result = coh_tmm('s', n_list, d_list, th_0, lam_vac)

    assert result['R'] == approx(np.array([0.06513963, 0.06122131]))


def test_coh_tmm_p_T():
    from solcore.absorption_calculator.tmm_core_vec import coh_tmm
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    result = coh_tmm('s', n_list, d_list, th_0, lam_vac)

    assert result['T'] == approx(np.array([1.15234466e-09, 4.13619185e-01]))


def test_coh_tmm_p_power_entering():
    from solcore.absorption_calculator.tmm_core_vec import coh_tmm
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    result = coh_tmm('p', n_list, d_list, th_0, lam_vac)

    assert result['power_entering'] == approx(np.array([0.96244989, 0.94994018]))


def test_coh_tmm_p_vw_list():
    from solcore.absorption_calculator.tmm_core_vec import coh_tmm
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    result = coh_tmm('p', n_list, d_list, th_0, lam_vac)

    assert result['vw_list'] == approx(np.array([[[ 0.00000000e+00+0.00000000e+00j,
          0.00000000e+00+0.00000000e+00j],
        [ 0.00000000e+00+0.00000000e+00j,
          0.00000000e+00+0.00000000e+00j]],
       [[ 1.17017431e+00-2.43748228e-01j,
          4.40679361e-02-1.53940000e-02j],
        [ 1.02922989e+00-7.82628087e-02j,
         -1.84573000e-01+1.79809491e-02j]],
       [[-8.59886075e-02+1.13689959e-01j,
          1.39851113e-09-3.01497601e-10j],
        [ 6.07730278e-01+5.07144030e-01j,
         -1.68609283e-01+8.64966880e-02j]],
       [[-1.23967610e-05+1.45623920e-05j,
          1.93199813e-06-2.24107827e-06j],
        [-6.52504472e-02+3.60299246e-01j,
         -3.33430797e-02-3.97852657e-02j]],
       [[ 1.82536479e-05-1.06422631e-05j,
          0.00000000e+00+0.00000000e+00j],
        [-9.02947159e-02-4.09448171e-01j,
          0.00000000e+00+0.00000000e+00j]]]))


def test_coh_tmm_kz_list():
    from solcore.absorption_calculator.tmm_core_vec import coh_tmm
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    result = coh_tmm('s', n_list, d_list, th_0, lam_vac)

    assert result['kz_list'] == approx(np.array([[0.02250959+0.00000000e+00j        , 0.00440866+0.00000000e+00j        ],
       [0.01435451+0.00687561j, 0.00404247+0.00074813j],
       [0.03118008+0.04748033j, 0.00515452+0.00110011j],
       [0.07823055+0.00000000e+00j      , 0.01413365+0.00000000e+00j       ],
       [0.06246792+0.01579948j, 0.01056188+0.00035793j]]), rel=1e-5)

def test_coh_tmm_th_list():
    from solcore.absorption_calculator.tmm_core_vec import coh_tmm
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    result = coh_tmm('s', n_list, d_list, th_0, lam_vac)

    assert result['th_list'] == approx(np.array([[0.3+0.j, 0.3+0.j],
       [0.38659626-0.16429512j, 0.3162772 -0.05459799j],
       [0.06789345-0.10235287j, 0.24849917-0.0507924j ],
       [0.08877261+0.j        , 0.09619234+0.j        ],
       [0.10445527-0.02621521j, 0.12826687-0.00429919j]]))


def test_coh_tmm_inputs():
    from solcore.absorption_calculator.tmm_core_vec import coh_tmm
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    result = coh_tmm('s', n_list, d_list, th_0, lam_vac)

    assert result['pol'] == 's'
    assert np.all(result['n_list'] == n_list)
    assert np.all(result['d_list'] == d_list[:,None])
    assert result['th_0'] == th_0
    assert np.all(result['lam_vac'] == lam_vac)

## end of tests for coh_tmm

def test_coh_tmm_reverse():
    from solcore.absorption_calculator.tmm_core_vec import coh_tmm_reverse
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    result = coh_tmm_reverse('s', n_list, d_list, th_0, lam_vac)

    assert result['r'] == approx(np.array([0.44495594-0.08125146j, 0.13483962-0.37537464j]))
    assert result['t'] == approx(np.array([ 5.64612420e-05-1.46509665e-05j, -1.94928443e-01-9.82812305e-01j]))
    assert result['R'] == approx(np.array([0.20458758, 0.15908784]))
    assert result['T'] == approx(np.array([1.22605928e-09, 4.19050975e-01]))
    assert result['power_entering'] == approx(np.array([0.75431194, 0.84091216]))
    assert result['vw_list'] == approx(np.array([[[ 0.00000000e+00+0.00000000e+00j,
          0.00000000e+00+0.00000000e+00j],
        [ 0.00000000e+00+0.00000000e+00j,
          0.00000000e+00+0.00000000e+00j]],
       [[ 9.35877457e-01+4.78628162e-02j,
          5.09078478e-01-1.29114277e-01j],
        [ 8.90678139e-01-4.74324242e-02j,
          2.44161483e-01-3.27942215e-01j]],
       [[-1.22416652e+00-1.12459390e-01j,
         -1.29313936e-08+6.59495568e-09j],
        [-1.13691217e+00+6.36142041e-01j,
         -2.07647974e-02-4.70097431e-02j]],
       [[-2.63661542e-04+6.60611182e-05j,
         -4.98027806e-07-4.91628647e-06j],
        [-1.06956388e+00-5.45729218e-01j,
          5.53546526e-02+6.76305892e-02j]],
       [[ 5.64612420e-05-1.46509665e-05j,
          0.00000000e+00+0.00000000e+00j],
        [-1.94928443e-01-9.82812305e-01j,
          0.00000000e+00+0.00000000e+00j]]]))
    assert result['kz_list'] == approx(np.array([[0.06246792+0.01579948j, 0.01056179+0.j        ],
       [0.07823055+0.j        , 0.01413365+0.j        ],
       [0.03118008+0.04748033j, 0.00515452+0.00110011j],
       [0.01435451+0.00687561j, 0.00404247+0.00074813j],
       [0.02250959+0.j        , 0.00440866+0.j        ]]), rel=1e-5)
    assert result['th_list'] == approx(np.array([[0.10445527-0.02621521j, 0.12841137+0.j        ],
       [0.08877261+0.j        , 0.09619234+0.j        ],
       [0.06789345-0.10235287j, 0.24849917-0.0507924j ],
       [0.38659626-0.16429512j, 0.3162772 -0.05459799j],
       [0.3       +0.j        , 0.3       +0.j        ]]))



def test_ellips_psi():
    from solcore.absorption_calculator.tmm_core_vec import ellips
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    assert ellips(n_list, d_list, th_0, lam_vac)['psi'] == approx(np.array([0.64939282, 0.73516374]))


def test_ellips_Delta():
    from solcore.absorption_calculator.tmm_core_vec import ellips

    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    assert ellips(n_list, d_list, th_0, lam_vac)['Delta'] == approx(np.array([0.09560384, 0.10886363]))


def test_unpolarized_RT_R():
    from solcore.absorption_calculator.tmm_core_vec import unpolarized_RT
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    assert unpolarized_RT(n_list, d_list, th_0, lam_vac)['R'] == approx(np.array([0.05134487, 0.05564057]))


def test_unpolarized_RT_T():
    from solcore.absorption_calculator.tmm_core_vec import unpolarized_RT
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    assert unpolarized_RT(n_list, d_list, th_0, lam_vac)['T'] == approx(np.array([1.19651603e-09, 4.17401823e-01]))


def test_position_resolved_s_poyn():
    from solcore.absorption_calculator.tmm_core_vec import position_resolved
    coh_tmm_data = {'r': np.array([0.14017645-0.2132843j , 0.22307786-0.10704008j]),
                    't': np.array([ 1.78669633e-05-9.79824244e-06j, -8.86075993e-02-4.05953564e-01j]),
                    'R': np.array([0.06513963, 0.06122131]),
                    'T': np.array([1.15234466e-09, 4.13619185e-01]),
                    'power_entering': np.array([0.93486037, 0.93877869]),
                    'vw_list': np.array([[[ 0.00000000e+00+0.00000000e+00j,
                                          0.00000000e+00+0.00000000e+00j],
                                        [ 0.00000000e+00+0.00000000e+00j,
                                          0.00000000e+00+0.00000000e+00j]],
                                       [[ 1.18358724e+00-2.33272105e-01j,
                                         -4.34107939e-02+1.99878010e-02j],
                                        [ 1.03160316e+00-7.28921467e-02j,
                                          1.91474694e-01-3.41479380e-02j]],
                                       [[-8.59535500e-02+1.06568462e-01j,
                                         -1.36521327e-09+2.83859953e-10j],
                                        [ 6.08369346e-01+5.06683493e-01j,
                                          1.75320349e-01-9.58306162e-02j]],
                                       [[-1.23112929e-05+1.37276841e-05j,
                                         -1.94390395e-06+2.16097082e-06j],
                                        [-6.54156818e-02+3.57104644e-01j,
                                          3.38453387e-02+4.04808706e-02j]],
                                       [[ 1.78669633e-05-9.79824244e-06j,
                                          0.00000000e+00+0.00000000e+00j],
                                        [-8.86075993e-02-4.05953564e-01j,
                                          0.00000000e+00+0.00000000e+00j]]]),
                    'kz_list':np.array([[0.02250959+0.j        , 0.00440866+0.j        ],
                                       [0.01435451+0.00687561j, 0.00404247+0.00074813j],
                                       [0.03118008+0.04748033j, 0.00515452+0.00110011j],
                                       [0.07823055+0.j        , 0.01413365+0.j        ],
                                       [0.06246792+0.01579948j, 0.01056188+0.00035793j]]),
                    'th_list':np.array([[0.3       +0.j        , 0.3       +0.j        ],
                                       [0.38659626-0.16429512j, 0.3162772 -0.05459799j],
                                       [0.06789345-0.10235287j, 0.24849917-0.0507924j ],
                                       [0.08877261+0.j        , 0.09619234+0.j        ],
                                       [0.10445527-0.02621521j, 0.12826687-0.00429919j]]),
                    'pol': 's',
                    'n_list':np.array([[1.5+0.j , 1.3+0.j ],
                                   [1. +0.4j, 1.2+0.2j],
                                   [2. +3.j , 1.5+0.3j],
                                   [5. +0.j , 4. +0.j ],
                                   [4. +1.j , 3. +0.1j]]),
                    'd_list':np.array([[np.inf], [ 200. ], [ 187.3], [1973.5], [np.inf]]),
                    'th_0': [0.3],
                    'lam_vac':np.array([400, 1770])}

    result = position_resolved([0, 1, 1, 2, 2, 3, 3, 4, 4], np.array([20, 10, 200, 20, 50.8, 10, 2000, 0, 200]), coh_tmm_data)

    assert result['poyn'] == approx(np.array([[0.00000000e+00, 8.24769450e-01, 2.59652919e-02, 3.88663596e-03,
        2.08618804e-04, 1.15234465e-09, 1.15234465e-09, 1.15234466e-09,
        2.07458676e-12],
       [0.00000000e+00, 9.18232659e-01, 6.12787705e-01, 5.74755463e-01,
        5.25103871e-01, 4.13619185e-01, 4.13619185e-01, 4.13619186e-01,
        3.58444455e-01]]))


def test_position_resolved_s_absor():
    from solcore.absorption_calculator.tmm_core_vec import position_resolved
    coh_tmm_data = {'r': np.array([0.14017645-0.2132843j , 0.22307786-0.10704008j]),
                    't': np.array([ 1.78669633e-05-9.79824244e-06j, -8.86075993e-02-4.05953564e-01j]),
                    'R': np.array([0.06513963, 0.06122131]),
                    'T': np.array([1.15234466e-09, 4.13619185e-01]),
                    'power_entering': np.array([0.93486037, 0.93877869]),
                    'vw_list': np.array([[[ 0.00000000e+00+0.00000000e+00j,
                                          0.00000000e+00+0.00000000e+00j],
                                        [ 0.00000000e+00+0.00000000e+00j,
                                          0.00000000e+00+0.00000000e+00j]],
                                       [[ 1.18358724e+00-2.33272105e-01j,
                                         -4.34107939e-02+1.99878010e-02j],
                                        [ 1.03160316e+00-7.28921467e-02j,
                                          1.91474694e-01-3.41479380e-02j]],
                                       [[-8.59535500e-02+1.06568462e-01j,
                                         -1.36521327e-09+2.83859953e-10j],
                                        [ 6.08369346e-01+5.06683493e-01j,
                                          1.75320349e-01-9.58306162e-02j]],
                                       [[-1.23112929e-05+1.37276841e-05j,
                                         -1.94390395e-06+2.16097082e-06j],
                                        [-6.54156818e-02+3.57104644e-01j,
                                          3.38453387e-02+4.04808706e-02j]],
                                       [[ 1.78669633e-05-9.79824244e-06j,
                                          0.00000000e+00+0.00000000e+00j],
                                        [-8.86075993e-02-4.05953564e-01j,
                                          0.00000000e+00+0.00000000e+00j]]]),
                    'kz_list':np.array([[0.02250959+0.j        , 0.00440866+0.j        ],
                                       [0.01435451+0.00687561j, 0.00404247+0.00074813j],
                                       [0.03118008+0.04748033j, 0.00515452+0.00110011j],
                                       [0.07823055+0.j        , 0.01413365+0.j        ],
                                       [0.06246792+0.01579948j, 0.01056188+0.00035793j]]),
                    'th_list':np.array([[0.3       +0.j        , 0.3       +0.j        ],
                                       [0.38659626-0.16429512j, 0.3162772 -0.05459799j],
                                       [0.06789345-0.10235287j, 0.24849917-0.0507924j ],
                                       [0.08877261+0.j        , 0.09619234+0.j        ],
                                       [0.10445527-0.02621521j, 0.12826687-0.00429919j]]),
                    'pol': 's',
                    'n_list':np.array([[1.5+0.j , 1.3+0.j ],
                                   [1. +0.4j, 1.2+0.2j],
                                   [2. +3.j , 1.5+0.3j],
                                   [5. +0.j , 4. +0.j ],
                                   [4. +1.j , 3. +0.1j]]),
                    'd_list':np.array([[np.inf], [ 200. ], [ 187.3], [1973.5], [np.inf]]),
                    'th_0': [0.3],
                    'lam_vac':np.array([400, 1770])}

    result = position_resolved([0, 1, 1, 2, 2, 3, 3, 4, 4], np.array([20, 10, 200, 20, 50.8, 10, 2000, 0, 200]),
                               coh_tmm_data)

    assert result['absor'] == approx(np.array([[0.00000000e+00, 1.02697381e-02, 1.64378595e-04, 3.69077590e-04,
        1.98105861e-05, 0.00000000e+00, 0.00000000e+00, 3.64128877e-11,
        6.55547748e-14],
       [0.00000000e+00, 2.04057536e-03, 1.07421979e-03, 1.78805499e-03,
        1.43715316e-03, 0.00000000e+00, 0.00000000e+00, 2.96091644e-04,
        2.56594499e-04]]))


def test_position_resolved_p_poyn():
    from solcore.absorption_calculator.tmm_core_vec import position_resolved
    coh_tmm_data = {'r':np.array([-0.12140058+0.15103645j, -0.21104259+0.07430242j]),
                    't':np.array([ 1.82536479e-05-1.06422631e-05j, -9.02947159e-02-4.09448171e-01j]),
                    'R':np.array([0.03755011, 0.05005982]),
                    'T':np.array([1.24068740e-09, 4.21184461e-01]),
                    'power_entering':np.array([0.96244989, 0.94994018]),
                    'vw_list':np.array([[[ 0.00000000e+00+0.00000000e+00j,
                                          0.00000000e+00+0.00000000e+00j],
                                        [ 0.00000000e+00+0.00000000e+00j,
                                          0.00000000e+00+0.00000000e+00j]],
                                       [[ 1.17017431e+00-2.43748228e-01j,
                                          4.40679361e-02-1.53940000e-02j],
                                        [ 1.02922989e+00-7.82628087e-02j,
                                         -1.84573000e-01+1.79809491e-02j]],
                                       [[-8.59886075e-02+1.13689959e-01j,
                                          1.39851113e-09-3.01497601e-10j],
                                        [ 6.07730278e-01+5.07144030e-01j,
                                         -1.68609283e-01+8.64966880e-02j]],
                                       [[-1.23967610e-05+1.45623920e-05j,
                                          1.93199813e-06-2.24107827e-06j],
                                        [-6.52504472e-02+3.60299246e-01j,
                                         -3.33430797e-02-3.97852657e-02j]],
                                       [[ 1.82536479e-05-1.06422631e-05j,
                                          0.00000000e+00+0.00000000e+00j],
                                        [-9.02947159e-02-4.09448171e-01j,
                                        0.00000000e+00+0.00000000e+00j]]]),
                    'kz_list':np.array([[0.02250959+0.j        , 0.00440866+0.j        ],
                                       [0.01435451+0.00687561j, 0.00404247+0.00074813j],
                                       [0.03118008+0.04748033j, 0.00515452+0.00110011j],
                                       [0.07823055+0.j        , 0.01413365+0.j        ],
                                       [0.06246792+0.01579948j, 0.01056188+0.00035793j]]),
                    'th_list':np.array([[0.3       +0.j        , 0.3       +0.j        ],
                                       [0.38659626-0.16429512j, 0.3162772 -0.05459799j],
                                       [0.06789345-0.10235287j, 0.24849917-0.0507924j ],
                                       [0.08877261+0.j        , 0.09619234+0.j        ],
                                       [0.10445527-0.02621521j, 0.12826687-0.00429919j]]),
                    'pol': 'p',
                    'n_list':np.array([[1.5+0.j , 1.3+0.j ],
                                       [1. +0.4j, 1.2+0.2j],
                                       [2. +3.j , 1.5+0.3j],
                                       [5. +0.j , 4. +0.j ],
                                       [4. +1.j , 3. +0.1j]]),
                    'd_list':np.array([[np.inf], [ 200. ], [ 187.3], [1973.5], [np.inf]]),
                    'th_0': [0.3],
                    'lam_vac':np.array([ 400, 1770])}

    result = position_resolved([0, 1, 1, 2, 2, 3, 3, 4, 4], np.array([20, 10, 200, 20, 50.8, 10, 2000, 0, 200]), coh_tmm_data)

    assert result['poyn'] == approx(np.array([[0.00000000e+00, 8.45520745e-01, 2.87382013e-02, 4.30170256e-03,
        2.30897897e-04, 1.24068740e-09, 1.24068740e-09, 1.24068740e-09,
        2.23363177e-12],
       [0.00000000e+00, 9.30639570e-01, 6.33536065e-01, 5.95958829e-01,
        5.45941168e-01, 4.21184460e-01, 4.21184460e-01, 4.21184461e-01,
        3.65000560e-01]]))


def test_position_resolved_p_absor():
    from solcore.absorption_calculator.tmm_core_vec import position_resolved
    coh_tmm_data = {'r':np.array([-0.12140058 + 0.15103645j, -0.21104259 + 0.07430242j]),
                    't':np.array([1.82536479e-05 - 1.06422631e-05j, -9.02947159e-02 - 4.09448171e-01j]),
                    'R':np.array([0.03755011, 0.05005982]),
                    'T':np.array([1.24068740e-09, 4.21184461e-01]),
                    'power_entering':np.array([0.96244989, 0.94994018]),
                    'vw_list':np.array([[[0.00000000e+00 + 0.00000000e+00j,
                                        0.00000000e+00 + 0.00000000e+00j],
                                       [0.00000000e+00 + 0.00000000e+00j,
                                        0.00000000e+00 + 0.00000000e+00j]],
                                      [[1.17017431e+00 - 2.43748228e-01j,
                                        4.40679361e-02 - 1.53940000e-02j],
                                       [1.02922989e+00 - 7.82628087e-02j,
                                        -1.84573000e-01 + 1.79809491e-02j]],
                                      [[-8.59886075e-02 + 1.13689959e-01j,
                                        1.39851113e-09 - 3.01497601e-10j],
                                       [6.07730278e-01 + 5.07144030e-01j,
                                        -1.68609283e-01 + 8.64966880e-02j]],
                                      [[-1.23967610e-05 + 1.45623920e-05j,
                                        1.93199813e-06 - 2.24107827e-06j],
                                       [-6.52504472e-02 + 3.60299246e-01j,
                                        -3.33430797e-02 - 3.97852657e-02j]],
                                      [[1.82536479e-05 - 1.06422631e-05j,
                                        0.00000000e+00 + 0.00000000e+00j],
                                       [-9.02947159e-02 - 4.09448171e-01j,
                                        0.00000000e+00 + 0.00000000e+00j]]]),
                    'kz_list':np.array([[0.02250959 + 0.j, 0.00440866 + 0.j],
                                      [0.01435451 + 0.00687561j, 0.00404247 + 0.00074813j],
                                      [0.03118008 + 0.04748033j, 0.00515452 + 0.00110011j],
                                      [0.07823055 + 0.j, 0.01413365 + 0.j],
                                      [0.06246792 + 0.01579948j, 0.01056188 + 0.00035793j]]),
                    'th_list':np.array([[0.3 + 0.j, 0.3 + 0.j],
                                      [0.38659626 - 0.16429512j, 0.3162772 - 0.05459799j],
                                      [0.06789345 - 0.10235287j, 0.24849917 - 0.0507924j],
                                      [0.08877261 + 0.j, 0.09619234 + 0.j],
                                      [0.10445527 - 0.02621521j, 0.12826687 - 0.00429919j]]),
                    'pol': 'p',
                    'n_list':np.array([[1.5 + 0.j, 1.3 + 0.j],
                                     [1. + 0.4j, 1.2 + 0.2j],
                                     [2. + 3.j, 1.5 + 0.3j],
                                     [5. + 0.j, 4. + 0.j],
                                     [4. + 1.j, 3. + 0.1j]]),
                    'd_list':np.array([[np.inf], [200.], [187.3], [1973.5], [np.inf]]),
                    'th_0': [0.3],
                    'lam_vac':np.array([400, 1770])}

    result = position_resolved([0, 1, 1, 2, 2, 3, 3, 4, 4], np.array([20, 10, 200, 20, 50.8, 10, 2000, 0, 200]),
                               coh_tmm_data)

    assert result['absor'] == approx(np.array([[0.00000000e+00, 1.08969642e-02, 5.17505408e-04, 4.08492602e-04,
        2.19262267e-05, 0.00000000e+00, 0.00000000e+00, 3.92044315e-11,
        7.05804410e-14],
       [0.00000000e+00, 1.91824256e-03, 1.12566565e-03, 1.77891531e-03,
        1.46979454e-03, 0.00000000e+00, 0.00000000e+00, 3.01509108e-04,
        2.61289301e-04]]))


def test_find_in_structure_exception():
    from solcore.absorption_calculator.tmm_core_vec import find_in_structure
    with raises(ValueError):
        find_in_structure([np.inf, 100, 200, np.inf], [0, 100, 200])


def test_find_in_structure():
    from solcore.absorption_calculator.tmm_core_vec import find_in_structure
    assert find_in_structure([200, 187.3,1973.5], np.linspace(0, 700, 10))[0] == approx(np.array([1, 1, 1, 2, 2, 3, 3, 3, 3, 3]))
    assert find_in_structure([200, 187.3,1973.5], np.linspace(0, 700, 10))[1] == approx(np.array([ 0., 77.77777778, 155.55555556, 33.33333333, 111.11111111,
                                                       1.58888889, 79.36666667,157.14444444, 234.92222222, 312.7]))


def test_find_in_structure_with_inf():
    from solcore.absorption_calculator.tmm_core_vec import find_in_structure_with_inf
    assert find_in_structure_with_inf([np.inf, 200, 187.3,1973.5, np.inf], np.linspace(0, 700, 10))[0] == approx(np.array([1, 1, 1, 2, 2, 3, 3, 3, 3, 3]))
    assert find_in_structure_with_inf([np.inf, 200, 187.3,1973.5, np.inf], np.linspace(0, 700, 10))[1] == approx(np.array([ 0., 77.77777778, 155.55555556, 33.33333333, 111.11111111,
                                                       1.58888889, 79.36666667,157.14444444, 234.92222222, 312.7]))



def test_layer_starts():
    from solcore.absorption_calculator.tmm_core_vec import layer_starts
    assert layer_starts([np.inf, 200, 187.3,1973.5, np.inf]) == approx(np.array([  -np.inf,    0. ,  200. ,  387.3, 2360.8]))


### tests for absorp_analytic_fn

def test_fill_in_s():
    from solcore.absorption_calculator.tmm_core_vec import absorp_analytic_fn
    coh_tmm_data = {'r': np.array([0.14017645-0.2132843j , 0.22307786-0.10704008j]),
                    't': np.array([ 1.78669633e-05-9.79824244e-06j, -8.86075993e-02-4.05953564e-01j]),
                    'R': np.array([0.06513963, 0.06122131]),
                    'T': np.array([1.15234466e-09, 4.13619185e-01]),
                    'power_entering': np.array([0.93486037, 0.93877869]),
                    'vw_list': np.array([[[ 0.00000000e+00+0.00000000e+00j,
                                          0.00000000e+00+0.00000000e+00j],
                                        [ 0.00000000e+00+0.00000000e+00j,
                                          0.00000000e+00+0.00000000e+00j]],
                                       [[ 1.18358724e+00-2.33272105e-01j,
                                         -4.34107939e-02+1.99878010e-02j],
                                        [ 1.03160316e+00-7.28921467e-02j,
                                          1.91474694e-01-3.41479380e-02j]],
                                       [[-8.59535500e-02+1.06568462e-01j,
                                         -1.36521327e-09+2.83859953e-10j],
                                        [ 6.08369346e-01+5.06683493e-01j,
                                          1.75320349e-01-9.58306162e-02j]],
                                       [[-1.23112929e-05+1.37276841e-05j,
                                         -1.94390395e-06+2.16097082e-06j],
                                        [-6.54156818e-02+3.57104644e-01j,
                                          3.38453387e-02+4.04808706e-02j]],
                                       [[ 1.78669633e-05-9.79824244e-06j,
                                          0.00000000e+00+0.00000000e+00j],
                                        [-8.86075993e-02-4.05953564e-01j,
                                          0.00000000e+00+0.00000000e+00j]]]),
                    'kz_list': np.array([[0.02250959+0.j        , 0.00440866+0.j        ],
                                       [0.01435451+0.00687561j, 0.00404247+0.00074813j],
                                       [0.03118008+0.04748033j, 0.00515452+0.00110011j],
                                       [0.07823055+0.j        , 0.01413365+0.j        ],
                                       [0.06246792+0.01579948j, 0.01056188+0.00035793j]]),
                    'th_list': np.array([[0.3       +0.j        , 0.3       +0.j        ],
                                       [0.38659626-0.16429512j, 0.3162772 -0.05459799j],
                                       [0.06789345-0.10235287j, 0.24849917-0.0507924j ],
                                       [0.08877261+0.j        , 0.09619234+0.j        ],
                                       [0.10445527-0.02621521j, 0.12826687-0.00429919j]]),
                    'pol': 's',
                    'n_list': np.array([[1.5+0.j , 1.3+0.j ],
                                   [1. +0.4j, 1.2+0.2j],
                                   [2. +3.j , 1.5+0.3j],
                                   [5. +0.j , 4. +0.j ],
                                   [4. +1.j , 3. +0.1j]]),
                    'd_list': np.array([[np.inf], [ 200. ], [ 187.3], [1973.5], [np.inf]]),
                    'th_0': [0.3],
                    'lam_vac': np.array([400, 1770])}

    a = absorp_analytic_fn().fill_in(coh_tmm_data, [1, 2])

    assert np.all([r == approx(e, rel=1e-5) for r, e in zip([a.a1, a.a3, a.A1, a.A2, a.A3],
                                                            [np.array([[0.01375122, 0.00149626],
                                                                    [0.09496066, 0.00220022]]),
                                                             np.array([[0.02870902, 0.00808494],
                                                                    [0.06236016, 0.01030904]]),
                                                             np.array([[2.00290348e-05, 5.19001441e-05],
                                                                    [2.55761667e-19, 1.02694503e-04]]),
                                                             np.array([[0.01276183, 0.00146736],
                                                                    [0.00246567, 0.00161252]]),
                                                             np.array([[-4.91455269e-04-1.18654706e-04j,  2.74416636e-04+2.91821819e-05j],
                                                            [ 1.94145098e-11-1.59280064e-11j,  1.49469559e-04+3.78492113e-04j]])]
                                           )])


def test_fill_in_p():
    from solcore.absorption_calculator.tmm_core_vec import absorp_analytic_fn
    coh_tmm_data = {'r': np.array([-0.12140058 + 0.15103645j, -0.21104259 + 0.07430242j]),
                    't': np.array([1.82536479e-05 - 1.06422631e-05j, -9.02947159e-02 - 4.09448171e-01j]),
                    'R': np.array([0.03755011, 0.05005982]),
                    'T': np.array([1.24068740e-09, 4.21184461e-01]),
                    'power_entering': np.array([0.96244989, 0.94994018]),
                    'vw_list': np.array([[[0.00000000e+00 + 0.00000000e+00j,
                                        0.00000000e+00 + 0.00000000e+00j],
                                       [0.00000000e+00 + 0.00000000e+00j,
                                        0.00000000e+00 + 0.00000000e+00j]],
                                      [[1.17017431e+00 - 2.43748228e-01j,
                                        4.40679361e-02 - 1.53940000e-02j],
                                       [1.02922989e+00 - 7.82628087e-02j,
                                        -1.84573000e-01 + 1.79809491e-02j]],
                                      [[-8.59886075e-02 + 1.13689959e-01j,
                                        1.39851113e-09 - 3.01497601e-10j],
                                       [6.07730278e-01 + 5.07144030e-01j,
                                        -1.68609283e-01 + 8.64966880e-02j]],
                                      [[-1.23967610e-05 + 1.45623920e-05j,
                                        1.93199813e-06 - 2.24107827e-06j],
                                       [-6.52504472e-02 + 3.60299246e-01j,
                                        -3.33430797e-02 - 3.97852657e-02j]],
                                      [[1.82536479e-05 - 1.06422631e-05j,
                                        0.00000000e+00 + 0.00000000e+00j],
                                       [-9.02947159e-02 - 4.09448171e-01j,
                                        0.00000000e+00 + 0.00000000e+00j]]]),
                    'kz_list': np.array([[0.02250959 + 0.j, 0.00440866 + 0.j],
                                      [0.01435451 + 0.00687561j, 0.00404247 + 0.00074813j],
                                      [0.03118008 + 0.04748033j, 0.00515452 + 0.00110011j],
                                      [0.07823055 + 0.j, 0.01413365 + 0.j],
                                      [0.06246792 + 0.01579948j, 0.01056188 + 0.00035793j]]),
                    'th_list': np.array([[0.3 + 0.j, 0.3 + 0.j],
                                      [0.38659626 - 0.16429512j, 0.3162772 - 0.05459799j],
                                      [0.06789345 - 0.10235287j, 0.24849917 - 0.0507924j],
                                      [0.08877261 + 0.j, 0.09619234 + 0.j],
                                      [0.10445527 - 0.02621521j, 0.12826687 - 0.00429919j]]),
                    'pol': 'p',
                    'n_list': np.array([[1.5 + 0.j, 1.3 + 0.j],
                                     [1. + 0.4j, 1.2 + 0.2j],
                                     [2. + 3.j, 1.5 + 0.3j],
                                     [5. + 0.j, 4. + 0.j],
                                     [4. + 1.j, 3. + 0.1j]]),
                    'd_list': np.array([[np.inf], [200.], [187.3], [1973.5], [np.inf]]),
                    'th_0': [0.3],
                    'lam_vac': np.array([400, 1770])}

    a = absorp_analytic_fn().fill_in(coh_tmm_data, [1, 2])

    assert np.all([r == approx(e, rel=1e-5) for r, e in zip([a.a1, a.a3, a.A1, a.A2, a.A3],
                                           [np.array([[0.01375122, 0.00149626],
                                                   [0.09496066, 0.00220022]]),
                                            np.array([[0.02870902, 0.00808494],
                                                    [0.06236016, 0.01030904]]),
                                            np.array([[2.01486783e-05, 4.74646908e-05],
                                                   [2.74885296e-19, 9.28559634e-05]]),
                                            np.array([[0.01321129, 0.00147049], [0.00272899, 0.00162005]]),
                                            np.array([[-3.47185512e-04 - 4.56403179e-05j, 2.11762306e-04 + 4.49397818e-06j],
                                                [2.01399935e-11 - 1.73429018e-11j, 1.32514817e-04 + 3.12222809e-04j]])]
                                           )])


def test_copy():
    from solcore.absorption_calculator.tmm_core_vec import absorp_analytic_fn
    a = absorp_analytic_fn()
    a.a1, a.a3, a.A1, a.A2, a.A3, a.d = 1, 0.5, 2, 7, 5+3*1j, [7, 3]
    b = a.copy()

    assert np.all([r == approx(e, rel=1e-5) for r, e in zip([a.a1, a.a3, a.A1, a.A2, a.A3],
                                                            [b.a1, b.a3, b.A1, b.A2, b.A3]
                                           )])


def test_run_array():
    from solcore.absorption_calculator.tmm_core_vec import absorp_analytic_fn
    a = absorp_analytic_fn()
    a.a1, a.a3, a.A1, a.A2, a.A3, a.d = np.array([0.01375122, 0.00149626]), \
                                        np.array([0.02870902, 0.00808494]),\
                                        np.array([2.00290348e-05, 5.19001441e-05]),\
                                        np.array([0.01276183, 0.00146736]),\
                                        np.array([-4.91455269e-04-1.18654706e-04j,  2.74416636e-04+2.91821819e-05j]), \
                                        np.array([200. ])

    assert a.run(np.linspace(0, 200, 7)) == approx(np.array([[0.01179895+0.j, 0.00772895+0.j, 0.00570666+0.j, 0.0043161 +0.j,
        0.00277534+0.j, 0.00118025+0.j, 0.00016438+0.j],
       [0.00206809+0.j, 0.00196401+0.j, 0.00182646+0.j, 0.00166052+0.j,
        0.00147356+0.j, 0.00127472+0.j, 0.00107422+0.j]]), rel=1e-4)


def test_run():
    from solcore.absorption_calculator.tmm_core_vec import absorp_analytic_fn
    a = absorp_analytic_fn()
    a.a1, a.a3, a.A1, a.A2, a.A3, a.d = np.array([0.01375122, 0.00149626]), \
                                        np.array([0.02870902, 0.00808494]),\
                                        np.array([2.00290348e-05, 5.19001441e-05]),\
                                        np.array([0.01276183, 0.00146736]),\
                                        np.array([-4.91455269e-04-1.18654706e-04j,  2.74416636e-04+2.91821819e-05j]), \
                                        np.array([200. ])

    assert a.run(200) == approx(np.array([0.00016438+0.j, 0.00107422+0.j]), rel=1e-4)


def test_scale():
    from solcore.absorption_calculator.tmm_core_vec import absorp_analytic_fn
    a = absorp_analytic_fn()
    a.a1, a.a3, a.A1, a.A2, a.A3, a.d = np.array([0.01375122, 0.00149626]), \
                                        np.array([0.02870902, 0.00808494]),\
                                        np.array([2.00290348e-05, 5.19001441e-05]),\
                                        np.array([0.01276183, 0.00146736]),\
                                        np.array([-4.91455269e-04-1.18654706e-04j,  2.74416636e-04+2.91821819e-05j]), \
                                        np.array([200. ])

    a.scale(0.7)
    assert np.all([r == approx(e, rel=1e-5) for r, e in zip([a.a1, a.a3, a.A1, a.A2, a.A3, a.d],
                                                     [np.array([0.01375122, 0.00149626]),np.array([0.02870902, 0.00808494]),
                                                      0.7*np.array([2.00290348e-05, 5.19001441e-05]), 0.7*np.array([0.01276183, 0.00146736]),
                                                      0.7*np.array([-4.91455269e-04-1.18654706e-04j,  2.74416636e-04+2.91821819e-05j]),
                                                      np.array([200.])]
                                                     )])



def test_add():
    from solcore.absorption_calculator.tmm_core_vec import absorp_analytic_fn
    a = absorp_analytic_fn()
    a.a1, a.a3, a.A1, a.A2, a.A3, a.d = np.array([0.01375122, 0.00149626]), \
                                        np.array([0.02870902, 0.00808494]),\
                                        np.array([2.00290348e-05, 5.19001441e-05]),\
                                        np.array([0.01276183, 0.00146736]),\
                                        np.array([-4.91455269e-04-1.18654706e-04j,  2.74416636e-04+2.91821819e-05j]), \
                                        np.array([200. ])

    b = absorp_analytic_fn()
    b.a1, b.a3, b.A1, b.A2, b.A3, b.d = np.array([0.01375122, 0.00149626]), \
                                        np.array([0.02870902, 0.00808494]),\
                                        np.array([2.e-05, 5e-05]),\
                                        np.array([0.05, 0.003]),\
                                        np.array([4.81455269e-04-1-04j,  3e-04+3-05j]), \
                                        np.array([200. ])

    a.add(b)

    assert np.all([r == approx(e, rel=1e-5) for r, e in zip([a.a1, a.a3, a.A1, a.A2, a.A3, a.d],
                                                     [np.array([0.01375122, 0.00149626]),
                                                      np.array([0.02870902, 0.00808494]),
                                                      np.array([2.00290348e-05, 5.19001441e-05]) + np.array([2.e-05, 5e-05]),
                                                      np.array([0.01276183, 0.00146736]) + np.array([0.05, 0.003]),
                                                      np.array([-4.91455269e-04 - 1.18654706e-04j,
                                                                      2.74416636e-04 + 2.91821819e-05j]) +  np.array([4.81455269e-04-1-04j,  3e-04+3-05j]),
                                                      np.array([200.])]
                                                     )])

def test_add_exception():
    from solcore.absorption_calculator.tmm_core_vec import absorp_analytic_fn
    a = absorp_analytic_fn()
    a.a1, a.a3, a.A1, a.A2, a.A3, a.d = np.array([0.01375122, 0.00149626]), \
                                        np.array([0.02870902, 0.00808494]),\
                                        np.array([2.00290348e-05, 5.19001441e-05]),\
                                        np.array([0.01276183, 0.00146736]),\
                                        np.array([-4.91455269e-04-1.18654706e-04j,  2.74416636e-04+2.91821819e-05j]), \
                                        np.array([200. ])

    b = absorp_analytic_fn()
    b.a1, b.a3, b.A1, b.A2, b.A3, b.d = np.array([7, 2]), \
                                        np.array([0.02870902, 0.00808494]),\
                                        np.array([2.e-05, 5e-05]),\
                                        np.array([0.05, 0.003]),\
                                        np.array([4.81455269e-04-1-04j,  3e-04+3-05j]), \
                                        np.array([200. ])

    with raises(ValueError):
        a.add(b)


def test_absorp_in_each_layer():
    from solcore.absorption_calculator.tmm_core_vec import absorp_in_each_layer
    coh_tmm_data = {'r': np.array([0.14017645-0.2132843j , 0.22307786-0.10704008j]),
                    't': np.array([ 1.78669633e-05-9.79824244e-06j, -8.86075993e-02-4.05953564e-01j]),
                    'R': np.array([0.06513963, 0.06122131]),
                    'T': np.array([1.15234466e-09, 4.13619185e-01]),
                    'power_entering': np.array([0.93486037, 0.93877869]),
                    'vw_list': np.array([[[ 0.00000000e+00+0.00000000e+00j,
                                          0.00000000e+00+0.00000000e+00j],
                                        [ 0.00000000e+00+0.00000000e+00j,
                                          0.00000000e+00+0.00000000e+00j]],
                                       [[ 1.18358724e+00-2.33272105e-01j,
                                         -4.34107939e-02+1.99878010e-02j],
                                        [ 1.03160316e+00-7.28921467e-02j,
                                          1.91474694e-01-3.41479380e-02j]],
                                       [[-8.59535500e-02+1.06568462e-01j,
                                         -1.36521327e-09+2.83859953e-10j],
                                        [ 6.08369346e-01+5.06683493e-01j,
                                          1.75320349e-01-9.58306162e-02j]],
                                       [[-1.23112929e-05+1.37276841e-05j,
                                         -1.94390395e-06+2.16097082e-06j],
                                        [-6.54156818e-02+3.57104644e-01j,
                                          3.38453387e-02+4.04808706e-02j]],
                                       [[ 1.78669633e-05-9.79824244e-06j,
                                          0.00000000e+00+0.00000000e+00j],
                                        [-8.86075993e-02-4.05953564e-01j,
                                          0.00000000e+00+0.00000000e+00j]]]),
                    'kz_list':np.array([[0.02250959+0.j        , 0.00440866+0.j        ],
                                       [0.01435451+0.00687561j, 0.00404247+0.00074813j],
                                       [0.03118008+0.04748033j, 0.00515452+0.00110011j],
                                       [0.07823055+0.j        , 0.01413365+0.j        ],
                                       [0.06246792+0.01579948j, 0.01056188+0.00035793j]]),
                    'th_list':np.array([[0.3       +0.j        , 0.3       +0.j        ],
                                       [0.38659626-0.16429512j, 0.3162772 -0.05459799j],
                                       [0.06789345-0.10235287j, 0.24849917-0.0507924j ],
                                       [0.08877261+0.j        , 0.09619234+0.j        ],
                                       [0.10445527-0.02621521j, 0.12826687-0.00429919j]]),
                    'pol': 's',
                    'n_list':np.array([[1.5+0.j , 1.3+0.j ],
                                   [1. +0.4j, 1.2+0.2j],
                                   [2. +3.j , 1.5+0.3j],
                                   [5. +0.j , 4. +0.j ],
                                   [4. +1.j , 3. +0.1j]]),
                    'd_list':np.array([[np.inf], [ 200. ], [ 187.3], [1973.5], [np.inf]]),
                    'th_0': [0.3],
                    'lam_vac':np.array([400, 1770])}

    assert absorp_in_each_layer(coh_tmm_data) == approx(np.array([[ 6.51396300e-02,  6.12213100e-02],
       [ 9.08895166e-01,  3.25991032e-01],
       [ 2.59652025e-02,  1.99168474e-01],
       [0, 0],
       [ 1.15234466e-09,  4.13619185e-01]]), abs=1e-10, rel=1e-6)


def test_inc_group_layers_exceptions():
    from solcore.absorption_calculator.tmm_core_vec import inc_group_layers
    n_list = np.array([[1.5, 1 + 0.4*1j, 2 + 3*1j, 5, 4+1*1j],
                      [1.3, 1.2 + 0.2*1j, 1.5 + 0.3*1j, 4, 3+0.1*1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    c_list = ['i', 'c', 'c', 'i', 'i']

    with raises(ValueError):
        inc_group_layers(n_list, np.array([d_list, d_list]), c_list)

    with raises(ValueError):
        inc_group_layers(n_list, d_list[1:-1], c_list)

    with raises(ValueError):
        inc_group_layers(n_list, d_list, ['c', 'c', 'c', 'i', 'i'])

    with raises(ValueError):
        inc_group_layers(n_list[1:-1], d_list, c_list)

    with raises(ValueError):
        inc_group_layers(n_list, d_list, ['i', 'c', 'c', 'f', 'i'])


def test_inc_group_layers():
    from solcore.absorption_calculator.tmm_core_vec import inc_group_layers
    n_list = np.array([[1.5, 1 + 0.4*1j, 2 + 3*1j, 5, 4+1*1j],
                      [1.3, 1.2 + 0.2*1j, 1.5 + 0.3*1j, 4, 3+0.1*1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    c_list = ['i', 'c', 'c', 'i', 'i']

    result = inc_group_layers(n_list, d_list, c_list)

    correct = [[np.inf, 200.0, 187.3, np.inf],
               np.array([[np.array([1.5+0.j, 1.3+0.j]), np.array([1. +0.4j, 1.2+0.2j]), np.array([2. +3.j , 1.5+0.3j]), np.array([5.+0.j, 4.+0.j])]]),
               [0, 3, 4], [0, np.nan, np.nan, 1, 2], [[0, 1, 2, 3]], [np.nan, [0, 1], [0, 2], np.nan, np.nan], [0],
               [np.nan, 0, np.nan], 1,  3,  5]

    assert result['stack_d_list'] == [[np.inf, 200.0, 187.3, np.inf]]
    assert result['stack_n_list'] == approx(np.array([[[1.5+0.j, 1.3+0.j], [1. +0.4j, 1.2+0.2j], [2. +3.j , 1.5+0.3j], [5.+0.j, 4.+0.j]]]))
    assert result['all_from_inc'] == [0, 3, 4]
    assert result['inc_from_all'] == [0, np.nan, np.nan, 1, 2]
    assert result['all_from_stack'] == [[0, 1, 2, 3]]
    assert result['stack_from_all'] == [np.nan, [0, 1], [0, 2], np.nan, np.nan]
    assert result['inc_from_stack'] == [0]
    assert result['stack_from_inc'] == [np.nan, 0, np.nan]
    assert result['num_stacks'] == 1
    assert result['num_inc_layers'] == 3
    assert result['num_layers'] == 5


## tests for inc_tmm

def test_inc_tmm_exception():
    from solcore.absorption_calculator.tmm_core_vec import inc_tmm

    with raises(ValueError):
        inc_tmm('s', [1+0.5*1j, 2, 3], [np.inf, 20, np.inf], ['i', 'i', 'i'], [0.5], 500)


def test_inc_tmm_s_R():
    from solcore.absorption_calculator.tmm_core_vec import inc_tmm
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    c_list = ['i', 'c', 'c', 'i', 'i']

    result = inc_tmm('s', n_list, d_list, c_list, th_0, lam_vac)

    assert result['R'] == approx(np.array([0.06513963, 0.09735299]))

def test_inc_tmm_s_R_incfirst():
    from solcore.absorption_calculator.tmm_core_vec import inc_tmm
    # testing the case where the coherent stack DOES NOT start right after the semi-infinite layer
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    c_list = ['i', 'i', 'c', 'i', 'i']

    result = inc_tmm('s', n_list, d_list, c_list, th_0, lam_vac)

    assert result['R'] == approx(np.array([0.08254069, 0.07335674]))


def test_inc_tmm_s_T():
    from solcore.absorption_calculator.tmm_core_vec import inc_tmm
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 1 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 100.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    c_list = ['i', 'c', 'i', 'i', 'i']

    result = inc_tmm('s', n_list, d_list, c_list, th_0, lam_vac)

    assert result['T'] == approx(np.array([1.22080402e-04, 3.88581656e-01]))


def test_inc_tmm_s_power_entering_list():
    from solcore.absorption_calculator.tmm_core_vec import inc_tmm
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 3 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 1973.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    c_list = ['i', 'c', 'c', 'i', 'i']

    result = inc_tmm('s', n_list, d_list, c_list, th_0, lam_vac)

    assert result['power_entering_list'] == approx(np.array([[1.00000000e+00, 1.09570873e-09, 1.09570873e-09],
       [1.00000000e+00, 3.87024124e-01, 3.87024124e-01]]))


def test_inc_tmm_s_VW_list():
    from solcore.absorption_calculator.tmm_core_vec import inc_tmm
    n_list = np.array([[1.5, 1 + 0.4 * 1j, 2 + 1 * 1j, 5, 4 + 1 * 1j],
                       [1.3, 1.2 + 0.2 * 1j, 1.5 + 0.3 * 1j, 4, 3 + 0.1 * 1j]]).T
    d_list = np.array([np.inf, 200, 187.3, 100.5, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([400, 1770])

    c_list = ['i', 'c', 'i', 'i', 'i']

    result = inc_tmm('s', n_list, d_list, c_list, th_0, lam_vac)

    assert result['VW_list'] == approx(np.array([[[np.nan,            np.nan],
        [4.99455825e-02, 6.91491754e-08],
        [1.25191080e-04, 3.11067761e-06],
        [1.22080402e-04, 0.00000000e+00]],
       [[np.nan,            np.nan],
        [7.30898359e-01, 7.46196606e-02],
        [3.96967199e-01, 8.38554284e-03],
        [3.88581656e-01, 0.00000000e+00]]]), nan_ok=True)



def test_inc_absorp_in_each_layer():
    from solcore.absorption_calculator.tmm_core_vec import inc_absorp_in_each_layer

    inc_data = {'T': np.array([1.22080402e-04, 3.88581656e-01]), 'R': np.array([0.0691011, 0.0934006]), 'VW_list': np.array([[[           np.nan,            np.nan],
        [4.99455825e-02, 6.91491754e-08],
        [1.25191080e-04, 3.11067761e-06],
        [1.22080402e-04, 0.00000000e+00]],
       [[           np.nan,            np.nan],
        [7.30898359e-01, 7.46196606e-02],
        [3.96967199e-01, 8.38554284e-03],
        [3.88581656e-01, 0.00000000e+00]]]), 'coh_tmm_data_list': [{'r': np.array([0.15815358-0.2099727j , 0.05062109-0.18394307j]), 't': np.array([-0.16908171+0.0889725j,  0.59979328+0.5149352j]), 'R': np.array([0.06910109, 0.03639755]), 'T': np.array([0.04994557, 0.73063361]), 'power_entering': np.array([0.93089891, 0.96360245]), 'vw_list': np.array([[[ 0.        +0.j        ,  0.        +0.j        ],
        [ 0.        +0.j        ,  0.        +0.j        ]],
       [[ 1.18009942-0.22823678j, -0.02194584+0.01826407j],
        [ 1.0438036 -0.08762503j,  0.00681749-0.09631804j]],
       [[-0.16908171+0.0889725j ,  0.        +0.j        ],
        [ 0.59979328+0.5149352j ,  0.        +0.j        ]]]), 'kz_list': np.array([[0.02250959+0.j        , 0.00440866+0.j        ],
       [0.01435451+0.00687561j, 0.00404247+0.00074813j],
       [0.03079749+0.01602339j, 0.00515452+0.00110011j]]), 'th_list': np.array([[0.3       +0.j        , 0.3       +0.j        ],
       [0.38659626-0.16429512j, 0.3162772 -0.05459799j],
       [0.17752825-0.08995035j, 0.24849917-0.0507924j ]]), 'pol': 's', 'n_list': np.array([[1.5+0.j , 1.3+0.j ],
       [1. +0.4j, 1.2+0.2j],
       [2. +1.j , 1.5+0.3j]]), 'd_list': np.array([[ np.inf],
       [200.],
       [ np.inf]]), 'th_0': np.array([0.3+0.j, 0.3+0.j]), 'lam_vac': np.array([ 400, 1770])}], 'coh_tmm_bdata_list': [{'r': np.array([0.36942774+0.02981787j, 0.05747204-0.01565208j]), 't': np.array([-0.29467151+0.00137131j,  0.57277308+0.75172209j]), 'R': np.array([0.13736596, 0.00354802]), 'T': np.array([0.06346552, 0.76391471]), 'power_entering': np.array([0.89366146, 0.98977083]), 'vw_list': np.array([[[ 0.        +0.j        ,  0.        +0.j        ],
        [ 0.        +0.j        ,  0.        +0.j        ]],
       [[ 1.37311562+0.0051288j , -0.00368788+0.02468907j],
        [ 1.13241606+0.01868049j, -0.07494402-0.03433257j]],
       [[-0.29467151+0.00137131j,  0.        +0.j        ],
        [ 0.57277308+0.75172209j,  0.        +0.j        ]]]), 'kz_list': np.array([[0.03079749+1.60233896e-02j, 0.00515452+1.10011324e-03j],
       [0.01435451+6.87561247e-03j, 0.00404247+7.48130214e-04j],
       [0.02250959+1.34865520e-19j, 0.00440866-1.52390418e-20j]]), 'th_list': np.array([[0.17752825-8.99503525e-02j, 0.24849917-5.07924048e-02j],
       [0.38659626-1.64295119e-01j, 0.3162772 -5.45979936e-02j],
       [0.3       -1.93687955e-17j, 0.3       +1.11743051e-17j]]), 'pol': 's', 'n_list': np.array([[2. +1.j , 1.5+0.3j],
       [1. +0.4j, 1.2+0.2j],
       [1.5+0.j , 1.3+0.j ]]), 'd_list': np.array([[ np.inf],
       [200.],
       [ np.inf]]), 'th_0': np.array([0.17752825-0.08995035j, 0.24849917-0.0507924j ]), 'lam_vac': np.array([ 400, 1770])}], 'stackFB_list': np.array([[[1.00000000e+00, 6.91491754e-08]],
       [[1.00000000e+00, 7.46196606e-02]]]), 'power_entering_list': np.array([[1.00000000e+00, 4.99455112e-02, 1.22080402e-04, 1.22080402e-04],
       [1.00000000e+00, 6.56777243e-01, 3.88581656e-01, 3.88581656e-01]]), 'stack_d_list': [[np.inf, 200.0, np.inf]], 'stack_n_list': [[np.array([1.5+0.j, 1.3+0.j]), np.array([1. +0.4j, 1.2+0.2j]), np.array([2. +1.j , 1.5+0.3j])]], 'all_from_inc': [0, 2, 3, 4], 'inc_from_all': [0, np.nan, 1, 2, 3], 'all_from_stack': [[0, 1, 2]], 'stack_from_all': [np.nan, [0, 1], np.nan, np.nan, np.nan], 'inc_from_stack': [0], 'stack_from_inc': [np.nan, 0, np.nan, np.nan], 'num_stacks': 1, 'num_inc_layers': 4, 'num_layers': 5}

    assert np.array(inc_absorp_in_each_layer(inc_data)) == approx(np.array([[6.91010944e-02, 9.34006064e-02],
       [8.80953397e-01, 2.49822147e-01],
       [4.98234308e-02, 2.68195587e-01],
       [0.00000000e+00, 0.00000000e+00],
       [1.22080402e-04, 3.88581656e-01]])
)


def test_inc_find_absorp_analytic_fn_exception():
    from solcore.absorption_calculator.tmm_core_vec import inc_find_absorp_analytic_fn

    inc_data = {'stack_from_all': [np.nan, [0, 1], np.nan, np.nan, np.nan]}

    with raises(ValueError):
        inc_find_absorp_analytic_fn(2, inc_data)


def test_inc_find_absorp_analytic_fn():
    from solcore.absorption_calculator.tmm_core_vec import inc_find_absorp_analytic_fn, inc_tmm

    n_list = np.array([[1, 2, 2, 1],
                      [1, 3, 3, 1]]).T
    d_list = np.array([np.inf, 100, 1000, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([100, 500])
    c_list = ['i', 'c', 'i', 'i']

    inc_tmm_data = inc_tmm('s', n_list, d_list, c_list, th_0, lam_vac)

    a = inc_find_absorp_analytic_fn(1, inc_tmm_data)

    assert np.all([r == approx(e, rel=1e-5) for r, e in zip([a.a1, a.a3, a.A1, a.A2, a.A3],
                                                            [np.array([0,0]),
                                                             np.array([0.24856865, 0.07503152]),
                                                             np.array([0, 0]),
                                                             np.array([0, 0]),
                                                             np.array([0, 0])])])



def test_inc_position_resolved():
    from solcore.absorption_calculator.tmm_core_vec import inc_position_resolved, inc_tmm, find_in_structure_with_inf

    n_list = np.array([[1, 2 + 0.5*1j, 2, 1],
                      [1, 3, 3 + 1*1j, 1]]).T
    d_list = np.array([np.inf, 100, 1000, np.inf])

    th_0 = [0.3]

    lam_vac = np.array([100, 500])
    c_list = ['i', 'c', 'i', 'i']

    inc_tmm_data = inc_tmm('s', n_list, d_list, c_list, th_0, lam_vac)

    dist = np.linspace(0, 1100, 12)

    layer, d_in_layer = find_in_structure_with_inf(d_list, dist)

    alphas = 4*np.pi*n_list.imag/lam_vac

    result = inc_position_resolved(layer, d_in_layer, inc_tmm_data, c_list, alphas)

    assert result == approx(np.array([[0.00000000e+00, 5.41619013e-02, 1.11248375e-04, 7.66705892e-03,
            4.38516795e+00, 2.50714861e+03, 1.43341858e+06, 8.19532120e+08,
            4.68553223e+11, 2.67887148e+14, 1.53159813e+17, 0.00000000e+00],
           [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
            0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
            0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.74764534e-13]]))



def test_beer_lambert():
    from solcore.absorption_calculator.tmm_core_vec import beer_lambert

    alphas = np.linspace(0, 1, 5)
    fraction = np.linspace(0.2, 1, 5)
    dist = np.linspace(0, 100e-9, 4)

    assert beer_lambert(alphas, fraction, dist)*1e9 == approx(np.array([[0, 0, 0, 0],
       [0.1, 0.1, 0.1, 0.1],
       [0.3, 0.3, 0.29999999, 0.29999999],
       [0.6, 0.59999999, 0.59999997, 0.59999996],
       [1, 0.99999997, 0.99999993, 0.9999999]]))
